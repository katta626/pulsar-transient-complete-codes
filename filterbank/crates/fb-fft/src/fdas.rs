use std::collections::HashMap;
use std::sync::Arc;

use fb_defs::{PeriodicCandidate, Series};
use rayon::prelude::*;
use rustfft::num_complex::{Complex, Complex32};
use rustfft::{Fft, FftPlanner};

use crate::{stats, FftSearchParams};

#[derive(Debug, Clone, Copy)]
pub struct FdasParams {
    pub zmax: f64,
    pub dz: f64,
    pub wmax: f64,
    pub dw: f64,
    pub refine: usize,
}

impl Default for FdasParams {
    fn default() -> Self {
        Self { zmax: 200.0, dz: 2.0, wmax: 0.0, dw: 20.0, refine: 50 }
    }
}

fn grid(max: f64, step: f64) -> Vec<f64> {
    if max <= 0.0 || step <= 0.0 {
        return vec![0.0];
    }
    let n = (max / step + 1e-9).floor() as i64;
    (-n..=n).map(|i| i as f64 * step).collect()
}

impl FdasParams {
    pub fn z_grid(&self) -> Vec<f64> {
        grid(self.zmax, self.dz)
    }

    pub fn w_grid(&self) -> Vec<f64> {
        grid(self.wmax, self.dw)
    }
}

pub fn halfwidth(z: f64, w: f64) -> usize {
    (z.abs() / 2.0 + w.abs() / 8.0).ceil() as usize + 16
}

pub fn response(z: f64, w: f64, shift: f64, hw: usize) -> Vec<Complex32> {
    let m = (8 * hw).max(512).next_power_of_two();
    let mut g: Vec<Complex<f64>> = (0..m)
        .map(|j| {
            let u = (j as f64 + 0.5) / m as f64;
            let s = u - 0.5;
            let ph = z / 2.0 * s * s + w / 6.0 * s * s * s - shift * u;
            let a = 2.0 * std::f64::consts::PI * ph.rem_euclid(1.0);
            Complex::new(a.cos(), a.sin())
        })
        .collect();
    FftPlanner::<f64>::new().plan_fft_forward(m).process(&mut g);
    let hw = hw as i64;
    (-hw..=hw)
        .map(|k| {
            let v = g[k.rem_euclid(m as i64) as usize];
            let a = -std::f64::consts::PI * k as f64 / m as f64;
            let r = v * Complex::new(a.cos(), a.sin()) / m as f64;
            Complex32::new(r.re as f32, r.im as f32)
        })
        .collect()
}

fn halfbin_kernel(z: f64, w: f64, hw: usize) -> (Vec<Complex32>, [f32; 2]) {
    let k0 = response(z, w, 0.0, hw);
    let k1 = response(z, w, 0.5, hw);
    let mut kh = vec![Complex32::new(0.0, 0.0); 4 * hw + 1];
    let mut norm = [0f32; 2];
    for i in 0..=2 * hw {
        kh[2 * i] = k0[i];
        norm[0] += k0[i].norm_sqr();
        if i < 2 * hw {
            kh[2 * i + 1] = k1[i];
            norm[1] += k1[i].norm_sqr();
        }
    }
    (kh, norm)
}

struct Engine {
    l: usize,
    hw: usize,
    u: usize,
    fwd: Arc<dyn Fft<f32>>,
    inv: Arc<dyn Fft<f32>>,
}

struct Kern {
    fft: Vec<Complex32>,
    norm: [f32; 2],
}

struct SubPlane {
    eng: Engine,
    map: Vec<usize>,
    kerns: Vec<(usize, Kern)>,
}

impl Engine {
    fn new(hw: usize) -> Self {
        let kl = 4 * hw + 1;
        let l = (4 * kl).next_power_of_two().max(2048);
        let mut p = FftPlanner::<f32>::new();
        Self { l, hw, u: l - kl + 1, fwd: p.plan_fft_forward(l), inv: p.plan_fft_inverse(l) }
    }

    fn with_output(hw: usize, u: usize) -> Self {
        let kl = 4 * hw + 1;
        let l = (u + kl - 1).next_power_of_two();
        let mut p = FftPlanner::<f32>::new();
        Self { l, hw, u: l - kl + 1, fwd: p.plan_fft_forward(l), inv: p.plan_fft_inverse(l) }
    }

    fn kernel(&self, z: f64, w: f64) -> Kern {
        let hw = halfwidth(z, w).min(self.hw);
        let (kh, norm) = halfbin_kernel(z, w, hw);
        let mut buf = vec![Complex32::new(0.0, 0.0); self.l];
        let off = 2 * (self.hw - hw);
        for (i, v) in kh.into_iter().enumerate() {
            buf[i + off] = v;
        }
        self.fwd.process(&mut buf);
        Kern { fft: buf, norm }
    }

    fn segment(&self, spec: &[Complex32], s0: i64) -> Vec<Complex32> {
        let mut buf: Vec<Complex32> = (0..self.l as i64)
            .map(|i| {
                let n = s0 + i;
                if n >= 0 && n % 2 == 0 && ((n / 2) as usize) < spec.len() {
                    spec[(n / 2) as usize]
                } else {
                    Complex32::new(0.0, 0.0)
                }
            })
            .collect();
        self.fwd.process(&mut buf);
        buf
    }

    fn correlate(&self, seg: &[Complex32], k: &Kern, s0: i64) -> Vec<f32> {
        let mut buf: Vec<Complex32> = seg.iter().zip(&k.fft).map(|(a, b)| a * b.conj()).collect();
        self.inv.process(&mut buf);
        let scale = 1.0 / (self.l as f32 * self.l as f32);
        let par = (s0.rem_euclid(2)) as usize;
        let g = [scale / k.norm[par].max(1e-30), scale / k.norm[1 - par].max(1e-30)];
        let mut out = vec![0f32; self.u];
        for (o, b) in out.chunks_mut(2).zip(buf.chunks(2)) {
            for (i, v) in o.iter_mut().enumerate() {
                *v = b[i].norm_sqr() * g[i];
            }
        }
        out
    }
}

fn stages(max_harm: u32) -> Vec<u32> {
    let mut v = vec![1u32];
    while v.last().unwrap() * 2 <= max_harm.max(1) {
        v.push(v.last().unwrap() * 2);
    }
    v
}

fn power_threshold(sigma: f64, numharm: u32, ntrials: f64) -> f32 {
    let mut p = numharm as f64;
    while stats::sigma_from_power(p, numharm, ntrials) < sigma && p < 1e7 {
        p *= 1.05;
    }
    p as f32
}

pub fn ntrials(nbins: usize, fp: &FdasParams) -> f64 {
    nbins as f64 * fp.z_grid().len() as f64 * fp.w_grid().len() as f64
}

fn search_range(nbins: usize, tsamp: f64, params: &FftSearchParams) -> (usize, usize) {
    let t = (nbins - 1) as f64 * 2.0 * tsamp;
    let hmax = stages(params.max_harm).last().copied().unwrap();
    let n_lo = ((2.0 * params.f_min * t).floor() as usize).max(2);
    let n_hi = ((2.0 * params.f_max * t * hmax as f64).ceil() as usize).min(2 * (nbins - 1));
    (n_lo, n_hi)
}

pub fn search(spec: &[Complex32], series: &Series, fp: &FdasParams, params: &FftSearchParams) -> Vec<PeriodicCandidate> {
    let nbins = spec.len();
    let t = (nbins - 1) as f64 * 2.0 * series.tsamp();
    let zs = fp.z_grid();
    let ws = fp.w_grid();
    let (nz, nw) = (zs.len(), ws.len());
    let eng = Engine::new(halfwidth(zs[nz - 1], ws[nw - 1]));
    let kerns: Vec<Kern> = (0..nz * nw).into_par_iter().map(|k| eng.kernel(zs[k / nw], ws[k % nw])).collect();
    let st = stages(params.max_harm);
    let (n_lo, n_hi) = search_range(nbins, series.tsamp(), params);
    let nt = ntrials((n_hi.saturating_sub(n_lo)) / 2, fp);
    let thr: Vec<f32> = st.iter().map(|&h| power_threshold(params.sigma_threshold, h, nt)).collect();
    let idx = |v: f64, step: f64, half: usize, len: usize| -> usize {
        if len == 1 {
            0
        } else {
            ((v / step).round() as i64 + half as i64).clamp(0, len as i64 - 1) as usize
        }
    };
    let mut subs: HashMap<(u32, u32), SubPlane> = HashMap::new();
    for &h in &st {
        for k in (1..h).step_by(2) {
            let q = k as f64 / h as f64;
            let map: Vec<usize> = (0..nz * nw)
                .map(|t| idx(zs[t / nw] * q, fp.dz, nz / 2, nz) * nw + idx(ws[t % nw] * q, fp.dw, nw / 2, nw))
                .collect();
            let mut need = map.clone();
            need.sort_unstable();
            need.dedup();
            let hw = need.iter().map(|&s| halfwidth(zs[s / nw], ws[s % nw])).max().unwrap();
            let e = Engine::with_output(hw, (eng.u * k as usize).div_ceil(h as usize) + 2);
            let kerns = need.par_iter().map(|&s| (s, e.kernel(zs[s / nw], ws[s % nw]))).collect();
            subs.insert((k, h), SubPlane { eng: e, map, kerns });
        }
    }
    let mid = series.header.start.add_seconds(t / 2.0);
    let obs_id = series.header.observation_id();
    let starts: Vec<usize> = (n_lo..n_hi).step_by(eng.u).collect();
    let bytes = 8 * nz * nw * eng.u;
    let batch = ((1usize << 30) / bytes).clamp(1, rayon::current_num_threads());
    let block = |n0: usize| -> Vec<PeriodicCandidate> {
        let mut out: Vec<PeriodicCandidate> = Vec::new();
        let ue = eng.u.min(n_hi - n0);
        let s0 = n0 as i64 - 2 * eng.hw as i64;
        let seg = eng.segment(spec, s0);
        let top: Vec<Vec<f32>> = kerns.par_iter().map(|k| eng.correlate(&seg, k, s0)).collect();
        let mut sum = top;
        for (si, &h) in st.iter().enumerate() {
            if h > 1 {
                for k in (1..h).step_by(2) {
                    let sub = &subs[&(k, h)];
                    let base = (n0 as u64 * k as u64 / h as u64) as i64;
                    let sq = base - 2 * sub.eng.hw as i64;
                    let segq = sub.eng.segment(spec, sq);
                    let planes: Vec<(usize, Vec<f32>)> = sub.kerns.par_iter().map(|(s, kn)| (*s, sub.eng.correlate(&segq, kn, sq))).collect();
                    let mut lookup: Vec<Option<&Vec<f32>>> = vec![None; nz * nw];
                    for (s, p) in &planes {
                        lookup[*s] = Some(p);
                    }
                    let src: Vec<u32> = (0..ue)
                        .map(|j| (((n0 + j) as u64 * k as u64 * 2 + h as u64) / (2 * h as u64) - base as u64) as u32)
                        .take_while(|&q| (q as usize) < sub.eng.u)
                        .collect();
                    sum.par_iter_mut().enumerate().for_each(|(t, row)| {
                        let plane = lookup[sub.map[t]].unwrap();
                        for (v, &q) in row.iter_mut().zip(&src) {
                            *v += plane[q as usize];
                        }
                    });
                }
            }
            let mut best = vec![(0f32, 0usize); ue];
            for (t, row) in sum.iter().enumerate() {
                for (b, &v) in best.iter_mut().zip(row) {
                    if v > b.0 {
                        *b = (v, t);
                    }
                }
            }
            let mut j = 0;
            while j < ue {
                if best[j].0 < thr[si] {
                    j += 1;
                    continue;
                }
                let mut b = j;
                while j < ue && best[j].0 >= thr[si] {
                    if best[j].0 > best[b].0 {
                        b = j;
                    }
                    j += 1;
                }
                let (pw, t_i) = best[b];
                let hf = h as f64;
                let f = (n0 + b) as f64 / 2.0 / t / hf;
                if f < params.f_min || f > params.f_max {
                    continue;
                }
                let (z, w) = (zs[t_i / nw] / hf, ws[t_i % nw] / hf);
                out.push(PeriodicCandidate {
                    obs_id: obs_id.clone(),
                    f_hz: f,
                    fdot: z / (t * t),
                    fddot: w / (t * t * t),
                    epoch: mid,
                    dm: series.dm,
                    snr: (2.0 * pw as f64).sqrt(),
                    sigma: stats::sigma_from_power(pw as f64, h, nt),
                    numharm: h,
                    method: format!("fdas z={z:.2} w={w:.2}"),
                    label: None,
                    probability: None,
                });
            }
        }
        out
    };
    let mut out: Vec<PeriodicCandidate> = Vec::new();
    for c in starts.chunks(batch) {
        out.extend(c.par_iter().flat_map_iter(|&n0| block(n0)).collect::<Vec<_>>());
    }
    out.sort_by(|a, b| b.sigma.partial_cmp(&a.sigma).unwrap());
    out.truncate(params.max_cands);
    out
}

pub fn power_at(spec: &[Complex32], r: f64, z: f64, w: f64) -> f64 {
    let hw = halfwidth(z, w);
    let m = r.floor();
    let k = response(z, w, m - r, hw);
    let mut acc = Complex::new(0f64, 0f64);
    let mut norm = 0f64;
    for (i, kv) in k.iter().enumerate() {
        let b = m as i64 + i as i64 - hw as i64;
        norm += kv.norm_sqr() as f64;
        if b > 0 && (b as usize) < spec.len() {
            let a = spec[b as usize];
            acc += Complex::new(a.re as f64, a.im as f64) * Complex::new(kv.re as f64, -kv.im as f64);
        }
    }
    acc.norm_sqr() / norm.max(1e-30)
}

pub fn harmonic_power(spec: &[Complex32], r: f64, z: f64, w: f64, numharm: u32) -> f64 {
    (1..=numharm).map(|h| power_at(spec, r * h as f64, z * h as f64, w * h as f64)).sum()
}

pub fn refine(spec: &[Complex32], c: &PeriodicCandidate, tsamp: f64, fp: &FdasParams, ntrials: f64) -> PeriodicCandidate {
    let t = (spec.len() - 1) as f64 * 2.0 * tsamp;
    let mut x = [c.f_hz * t, c.fdot * t * t, c.fddot * t * t * t];
    let use_w = fp.wmax > 0.0;
    let mut step = [0.5, fp.dz.max(0.5) / 2.0, if use_w { fp.dw.max(1.0) / 2.0 } else { 0.0 }];
    let h = c.numharm.max(1);
    let f = |p: &[f64; 3]| harmonic_power(spec, p[0], p[1], p[2], h);
    let mut best = f(&x);
    for _ in 0..200 {
        let mut moved = false;
        for d in 0..3 {
            if step[d] == 0.0 {
                continue;
            }
            for sgn in [1.0, -1.0] {
                let mut y = x;
                y[d] += sgn * step[d];
                let v = f(&y);
                if v > best {
                    best = v;
                    x = y;
                    moved = true;
                    break;
                }
            }
        }
        if !moved {
            for s in step.iter_mut() {
                *s /= 2.0;
            }
            if step[0] < 1.0 / 128.0 {
                break;
            }
        }
    }
    let mut o = c.clone();
    o.f_hz = x[0] / t;
    o.fdot = x[1] / (t * t);
    o.fddot = x[2] / (t * t * t);
    o.snr = (2.0 * best).sqrt();
    o.sigma = stats::sigma_from_power(best, h, ntrials);
    o.method = format!("fdas z={:.2} w={:.2} refined", x[1], x[2]);
    o
}

pub fn search_series(series: &Series, birdies: &[(f64, f64)], fp: &FdasParams, params: &FftSearchParams) -> Vec<PeriodicCandidate> {
    let spec = crate::prepare_width(series, birdies, 100.max((8.0 * fp.zmax) as usize));
    let mut c = search(&spec, series, fp, params);
    let n = fp.refine.min(c.len());
    let (n_lo, n_hi) = search_range(spec.len(), series.tsamp(), params);
    let nt = ntrials(n_hi.saturating_sub(n_lo) / 2, fp);
    let refined: Vec<PeriodicCandidate> = c[..n].par_iter().map(|x| refine(&spec, x, series.tsamp(), fp, nt)).collect();
    c.splice(..n, refined);
    c.retain(|x| x.sigma >= params.sigma_threshold);
    c.sort_by(|a, b| b.sigma.partial_cmp(&a.sigma).unwrap());
    c
}
