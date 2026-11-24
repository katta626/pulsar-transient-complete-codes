pub mod stats;

use fb_defs::{Epoch, PeriodicCandidate, Series, SPEED_OF_LIGHT};
use rayon::prelude::*;
use realfft::RealFftPlanner;
use rustfft::num_complex::Complex32;

pub fn good_size(n: usize) -> usize {
    let mut m = n.max(2);
    loop {
        let mut k = m;
        for p in [2usize, 3, 5, 7] {
            while k % p == 0 {
                k /= p;
            }
        }
        if k == 1 {
            return m;
        }
        m += 1;
    }
}

pub fn rfft(x: &[f32]) -> Vec<Complex32> {
    let n = x.len();
    let mut planner = RealFftPlanner::<f32>::new();
    let r = planner.plan_fft_forward(n);
    let mut inp = x.to_vec();
    let mut out = r.make_output_vec();
    r.process(&mut inp, &mut out).expect("fft");
    out
}

pub fn presto_pack(spec: &[Complex32]) -> Vec<Complex32> {
    let n = spec.len() - 1;
    let mut out = spec[..n].to_vec();
    out[0] = Complex32::new(spec[0].re, spec[n].re);
    out
}

pub fn powers(spec: &[Complex32]) -> Vec<f32> {
    spec.iter().map(|c| c.norm_sqr()).collect()
}

pub fn deredden(spec: &mut [Complex32], start_width: usize, end_width: usize) {
    let n = spec.len();
    let mut lo = 1usize;
    let mut w = start_width.max(6);
    let mut meds: Vec<(f64, f32)> = Vec::new();
    while lo < n {
        let hi = (lo + w).min(n);
        let mut p: Vec<f32> = spec[lo..hi].iter().map(|c| c.norm_sqr()).collect();
        let mid = p.len() / 2;
        p.select_nth_unstable_by(mid, |a, b| a.partial_cmp(b).unwrap());
        let med = (p[mid] / std::f32::consts::LN_2).max(1e-30);
        meds.push((0.5 * (lo + hi) as f64, med));
        lo = hi;
        w = ((w as f64 * (1.0 + (w as f64).ln() / 6.0)) as usize).max(w + 1).min(end_width.max(start_width));
    }
    let mut k = 0usize;
    for (i, c) in spec.iter_mut().enumerate().skip(1) {
        let x = i as f64;
        while k + 1 < meds.len() && meds[k + 1].0 < x {
            k += 1;
        }
        let m = if k + 1 < meds.len() && x > meds[k].0 {
            let (x0, m0) = meds[k];
            let (x1, m1) = meds[k + 1];
            let f = ((x - x0) / (x1 - x0)) as f32;
            m0 * (1.0 - f) + m1 * f
        } else {
            meds[k].1
        };
        *c /= m.sqrt();
    }
    spec[0] = Complex32::new(0.0, 0.0);
}

pub fn interbin(spec: &[Complex32]) -> Vec<f32> {
    let n = spec.len();
    let mut out = vec![0f32; 2 * n];
    let k = std::f32::consts::PI / 4.0;
    for i in 0..n {
        out[2 * i] = spec[i].norm_sqr();
        if i + 1 < n {
            out[2 * i + 1] = ((spec[i] - spec[i + 1]) * k).norm_sqr();
        }
    }
    out
}

pub fn zap(spec: &mut [Complex32], t_obs: f64, birdies: &[(f64, f64)]) {
    for &(f, w) in birdies {
        let lo = (((f - w / 2.0) * t_obs).floor().max(0.0)) as usize;
        let hi = (((f + w / 2.0) * t_obs).ceil() as usize).min(spec.len());
        for c in &mut spec[lo..hi] {
            *c = Complex32::new(0.0, 0.0);
        }
    }
}

pub fn harmonic_sums_exact(p: &[f32], max_harm: u32) -> Vec<(u32, Vec<f32>)> {
    let n = p.len();
    let mut out = Vec::new();
    let mut h = 1u32;
    while h <= max_harm {
        let mut s = vec![0f32; n];
        s.par_iter_mut().enumerate().for_each(|(i, v)| {
            if i == 0 {
                return;
            }
            let mut acc = 0f32;
            for k in 1..=h {
                let src = ((i as f64) * k as f64 / h as f64).round() as usize;
                if src < n {
                    acc += p[src];
                }
            }
            *v = acc;
        });
        out.push((h, s));
        h *= 2;
    }
    out
}

#[derive(Debug, Clone, Copy)]
pub struct FftSearchParams {
    pub max_harm: u32,
    pub sigma_threshold: f64,
    pub f_min: f64,
    pub f_max: f64,
    pub max_cands: usize,
}

impl Default for FftSearchParams {
    fn default() -> Self {
        Self { max_harm: 16, sigma_threshold: 6.0, f_min: 0.1, f_max: 2000.0, max_cands: 200 }
    }
}

fn peaks(sum: &[f32], numharm: u32, bins_per_hz: f64, params: &FftSearchParams, ntrials: f64, half_bins: bool) -> Vec<(usize, f64, f32)> {
    let fac = if half_bins { 2.0 } else { 1.0 } * numharm as f64;
    let lo = ((params.f_min * bins_per_hz * fac) as usize).max(1);
    let hi = ((params.f_max * bins_per_hz * fac) as usize).min(sum.len().saturating_sub(1));
    let mut out = Vec::new();
    if lo >= hi {
        return out;
    }
    let thr_power = {
        let mut p = numharm as f64;
        while stats::sigma_from_power(p, numharm, ntrials) < params.sigma_threshold && p < 1e6 {
            p *= 1.1;
        }
        p as f32
    };
    let mut i = lo;
    while i < hi {
        if sum[i] >= thr_power {
            let mut j = i;
            let mut best = i;
            while j < hi && sum[j] >= thr_power {
                if sum[j] > sum[best] {
                    best = j;
                }
                j += 1;
            }
            let sig = stats::sigma_from_power(sum[best] as f64, numharm, ntrials);
            out.push((best, sig, sum[best]));
            i = j;
        } else {
            i += 1;
        }
    }
    out
}

pub fn search_spectrum(spec: &[Complex32], series: &Series, accel: f64, params: &FftSearchParams) -> Vec<PeriodicCandidate> {
    let n_t = (spec.len() - 1) * 2;
    let t_obs = n_t as f64 * series.tsamp();
    let ib = interbin(spec);
    let sums = harmonic_sums_exact(&ib, params.max_harm);
    let ntrials = spec.len() as f64 * 2.0;
    let mid = series.header.start.add_seconds(t_obs / 2.0);
    let mut out = Vec::new();
    for (h, s) in &sums {
        for (bin, sig, pw) in peaks(s, *h, t_obs, params, ntrials, true) {
            let f = bin as f64 / 2.0 / t_obs / *h as f64;
            let fdot = -accel * f / SPEED_OF_LIGHT;
            out.push(PeriodicCandidate {
                obs_id: series.header.observation_id(),
                f_hz: f,
                fdot,
                fddot: 0.0,
                epoch: mid,
                dm: series.dm,
                snr: (2.0 * pw as f64).sqrt(),
                sigma: sig,
                numharm: *h,
                method: if accel == 0.0 { "fft".into() } else { format!("fft accel={accel:.3}") },
                label: None,
                probability: None,
            });
        }
    }
    out.sort_by(|a, b| b.sigma.partial_cmp(&a.sigma).unwrap());
    out.truncate(params.max_cands);
    out
}

pub fn prepare(series: &Series, birdies: &[(f64, f64)]) -> Vec<Complex32> {
    let n = good_size(series.len());
    let (mean, _) = series.mean_std();
    let mut x: Vec<f32> = series.data.iter().map(|v| *v - mean as f32).collect();
    x.resize(n, 0.0);
    let mut spec = rfft(&x);
    let t_obs = n as f64 * series.tsamp();
    zap(&mut spec, t_obs, birdies);
    deredden(&mut spec, 6, 100);
    spec
}

pub fn resample(x: &[f32], tsamp: f64, accel: f64) -> Vec<f32> {
    let n = x.len();
    let k = accel * tsamp / (2.0 * SPEED_OF_LIGHT);
    (0..n)
        .map(|i| {
            let fi = i as f64;
            let j = (fi + k * fi * (fi - n as f64)).round();
            if j >= 0.0 && (j as usize) < n { x[j as usize] } else { 0.0 }
        })
        .collect()
}

pub fn accel_trials(t_obs: f64, tsamp: f64, max_accel: f64, pulse_width: f64, tol: f64) -> Vec<f64> {
    if max_accel <= 0.0 {
        return vec![0.0];
    }
    let tw = pulse_width.max(tsamp);
    let step = 2.0 * tw * 24.0 * SPEED_OF_LIGHT / (t_obs * t_obs) * (tol * tol - 1.0).sqrt();
    let mut out = vec![0.0];
    let mut a = step;
    while a <= max_accel {
        out.push(a);
        out.push(-a);
        a += step;
    }
    out
}

pub fn search_series(series: &Series, accels: &[f64], birdies: &[(f64, f64)], params: &FftSearchParams) -> Vec<PeriodicCandidate> {
    let base = prepare(series, birdies);
    let mut all: Vec<PeriodicCandidate> = accels
        .par_iter()
        .flat_map(|&a| {
            if a == 0.0 {
                search_spectrum(&base, series, 0.0, params)
            } else {
                let n = good_size(series.len());
                let (mean, _) = series.mean_std();
                let mut x: Vec<f32> = series.data.iter().map(|v| *v - mean as f32).collect();
                x.resize(n, 0.0);
                let r = resample(&x, series.tsamp(), a);
                let mut spec = rfft(&r);
                zap(&mut spec, n as f64 * series.tsamp(), birdies);
                deredden(&mut spec, 6, 100);
                search_spectrum(&spec, series, a, params)
            }
        })
        .collect();
    all.sort_by(|a, b| b.sigma.partial_cmp(&a.sigma).unwrap());
    all
}

pub fn is_harmonic(f1: f64, f2: f64, rel_tol: f64, max_ratio: u32) -> bool {
    let (lo, hi) = if f1 < f2 { (f1, f2) } else { (f2, f1) };
    for den in 1..=max_ratio {
        let top = if den == 1 { 64.max(max_ratio) } else { max_ratio * den };
        for num in den..=top {
            let r = num as f64 / den as f64;
            if (hi / (lo * r) - 1.0).abs() < rel_tol {
                return true;
            }
        }
    }
    false
}

pub fn sift(cands: &[PeriodicCandidate], t_obs: f64, min_dm_hits: usize, max_ratio: u32) -> Vec<PeriodicCandidate> {
    let mut c = cands.to_vec();
    c.sort_by(|a, b| b.sigma.partial_cmp(&a.sigma).unwrap());
    let df = 1.5 / t_obs;
    let mut groups: Vec<(PeriodicCandidate, Vec<f64>)> = Vec::new();
    for x in c {
        let mut placed = false;
        for g in groups.iter_mut() {
            if (g.0.f_hz - x.f_hz).abs() < df {
                if !g.1.iter().any(|d| (*d - x.dm).abs() < 1e-9) {
                    g.1.push(x.dm);
                }
                placed = true;
                break;
            }
        }
        if !placed {
            groups.push((x.clone(), vec![x.dm]));
        }
    }
    let mut kept: Vec<PeriodicCandidate> = Vec::new();
    for (g, dms) in groups {
        if dms.len() < min_dm_hits {
            continue;
        }
        if kept.iter().any(|k| is_harmonic(k.f_hz, g.f_hz, df / g.f_hz.max(1e-9) + 2e-4, max_ratio)) {
            continue;
        }
        kept.push(g);
    }
    kept
}

pub fn epoch_mid(start: &Epoch, n: usize, tsamp: f64) -> Epoch {
    start.add_seconds(n as f64 * tsamp / 2.0)
}
