use fb_defs::{PeriodicCandidate, Series};
use rayon::prelude::*;

pub fn transform(x: &[f32], m: usize, p: usize) -> Vec<f32> {
    if m == 1 {
        return x[..p].to_vec();
    }
    let h = m / 2;
    let t = m - h;
    let head = transform(&x[..h * p], h, p);
    let tail = transform(&x[h * p..m * p], t, p);
    let mut out = vec![0f32; m * p];
    let denom = (m - 1) as f64;
    for s in 0..m {
        let kh = if h > 1 { ((s as f64) * (h - 1) as f64 / denom).round() as usize } else { 0 };
        let kt = if t > 1 { ((s as f64) * (t - 1) as f64 / denom).round() as usize } else { 0 };
        let shift = ((s as f64) * h as f64 / denom).round() as usize % p;
        let hr = &head[kh * p..(kh + 1) * p];
        let tr = &tail[kt * p..(kt + 1) * p];
        let o = &mut out[s * p..(s + 1) * p];
        for j in 0..p {
            o[j] = hr[j] + tr[(j + shift) % p];
        }
    }
    out
}

pub fn boxcar_snr(profile: &[f32], widths: &[usize], sigma: f64) -> (f64, usize) {
    let n = profile.len();
    let total: f64 = profile.iter().map(|v| *v as f64).sum();
    let mut best = (f64::NEG_INFINITY, 1);
    for &w in widths {
        if w >= n {
            continue;
        }
        let mut s: f64 = profile[..w].iter().map(|v| *v as f64).sum();
        let mut mx = s;
        for i in 0..n {
            s += profile[(i + w) % n] as f64 - profile[i] as f64;
            mx = mx.max(s);
        }
        let nf = n as f64;
        let wf = w as f64;
        let hgt = ((nf - wf) / (nf * wf)).sqrt();
        let b = wf * hgt / (nf - wf);
        let snr = ((hgt + b) * mx - b * total) / sigma;
        if snr > best.0 {
            best = (snr, w);
        }
    }
    best
}

pub fn downsample(x: &[f32], factor: f64) -> Vec<f32> {
    if factor <= 1.0 {
        return x.to_vec();
    }
    let n = (x.len() as f64 / factor).floor() as usize;
    (0..n)
        .map(|i| {
            let a = i as f64 * factor;
            let b = a + factor;
            let mut s = 0f64;
            let mut k = a;
            while k < b - 1e-9 {
                let idx = k.floor() as usize;
                let next = (idx as f64 + 1.0).min(b);
                s += x.get(idx).copied().unwrap_or(0.0) as f64 * (next - k);
                k = next;
            }
            (s / factor.sqrt()) as f32
        })
        .collect()
}

#[derive(Debug, Clone, Copy)]
pub struct FfaParams {
    pub period_min: f64,
    pub period_max: f64,
    pub bins_min: usize,
    pub bins_max: usize,
    pub snr_threshold: f64,
    pub duty_max: f64,
}

impl Default for FfaParams {
    fn default() -> Self {
        Self { period_min: 0.5, period_max: 30.0, bins_min: 240, bins_max: 260, snr_threshold: 7.0, duty_max: 0.2 }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct FfaPoint {
    pub period: f64,
    pub snr: f64,
    pub width: usize,
    pub bins: usize,
}

fn detrend_normalise(x: &[f32], window: usize) -> Vec<f32> {
    let n = x.len();
    let w = window.max(16).min(n.max(1));
    let nb = n.div_ceil(w);
    let mut out = x.to_vec();
    let meds: Vec<f32> = (0..nb)
        .map(|b| {
            let mut v = x[b * w..((b + 1) * w).min(n)].to_vec();
            let mid = v.len() / 2;
            v.select_nth_unstable_by(mid, |a, b| a.partial_cmp(b).unwrap());
            v[mid]
        })
        .collect();
    for i in 0..n {
        out[i] -= meds[(i / w).min(nb - 1)];
    }
    let var = out.iter().map(|v| (*v as f64).powi(2)).sum::<f64>() / n.max(1) as f64;
    let sd = var.sqrt().max(1e-12) as f32;
    for v in out.iter_mut() {
        *v /= sd;
    }
    out
}

pub fn periodogram(series: &Series, params: &FfaParams) -> Vec<FfaPoint> {
    let tsamp = series.tsamp();
    let x = detrend_normalise(&series.data, ((params.period_max * 3.0) / tsamp) as usize);
    let mut stages: Vec<f64> = Vec::new();
    let mut ds = params.period_min / (params.bins_min as f64 * tsamp);
    let geo = params.bins_max as f64 / params.bins_min as f64;
    while ds * tsamp * params.bins_min as f64 <= params.period_max {
        stages.push(ds.max(1.0));
        ds *= geo;
    }
    stages
        .par_iter()
        .flat_map(|&ds| {
            let y = downsample(&x, ds);
            let dt = ds * tsamp;
            let mut pts = Vec::new();
            for b in params.bins_min..params.bins_max {
                let m = y.len() / b;
                if m < 2 {
                    continue;
                }
                let maxw = ((b as f64 * params.duty_max) as usize).max(1);
                let mut widths = vec![1usize];
                while *widths.last().unwrap() < maxw {
                    let l = *widths.last().unwrap();
                    widths.push(((l as f64 * 1.5).ceil() as usize).max(l + 1));
                }
                let t = transform(&y[..m * b], m, b);
                let sigma = (m as f64).sqrt();
                for s in 0..m {
                    let prof = &t[s * b..(s + 1) * b];
                    let (snr, w) = boxcar_snr(prof, &widths, sigma);
                    let period = (b as f64 + s as f64 / (m - 1).max(1) as f64) * dt;
                    if period >= params.period_min && period <= params.period_max {
                        pts.push(FfaPoint { period, snr, width: w, bins: b });
                    }
                }
            }
            pts
        })
        .collect()
}

pub fn search(series: &Series, params: &FfaParams) -> Vec<PeriodicCandidate> {
    let mut pts = periodogram(series, params);
    pts.sort_by(|a, b| a.period.partial_cmp(&b.period).unwrap());
    let t_obs = series.len() as f64 * series.tsamp();
    let mut out: Vec<PeriodicCandidate> = Vec::new();
    let mut i = 0;
    while i < pts.len() {
        if pts[i].snr < params.snr_threshold {
            i += 1;
            continue;
        }
        let mut best = pts[i];
        let mut j = i;
        while j < pts.len() && pts[j].snr >= params.snr_threshold && (pts[j].period - best.period).abs() < best.period * best.period / t_obs {
            if pts[j].snr > best.snr {
                best = pts[j];
            }
            j += 1;
        }
        out.push(PeriodicCandidate {
            obs_id: series.header.observation_id(),
            f_hz: 1.0 / best.period,
            fdot: 0.0,
            fddot: 0.0,
            epoch: series.header.start,
            dm: series.dm,
            snr: best.snr,
            sigma: best.snr,
            numharm: 0,
            method: format!("ffa width={}/{}", best.width, best.bins),
            label: None,
            probability: None,
        });
        i = j.max(i + 1);
    }
    out.sort_by(|a, b| b.snr.partial_cmp(&a.snr).unwrap());
    out
}
