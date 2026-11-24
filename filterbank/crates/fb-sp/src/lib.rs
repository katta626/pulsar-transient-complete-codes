use fb_defs::{Epoch, ObsHeader, SinglePulseCandidate};
use rayon::prelude::*;

fn median(v: &mut [f32]) -> f32 {
    if v.is_empty() {
        return 0.0;
    }
    let mid = v.len() / 2;
    v.select_nth_unstable_by(mid, |a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    v[mid]
}

pub fn baseline_subtract(x: &mut [f32], window: usize) {
    let n = x.len();
    let w = window.max(8).min(n.max(1));
    let nb = n.div_ceil(w);
    let meds: Vec<f32> = (0..nb).map(|b| median(&mut x[b * w..((b + 1) * w).min(n)].to_vec())).collect();
    for i in 0..n {
        let pos = (i as f64 + 0.5) / w as f64 - 0.5;
        let b0 = pos.floor().clamp(0.0, (nb - 1) as f64) as usize;
        let b1 = (b0 + 1).min(nb - 1);
        let f = (pos - b0 as f64).clamp(0.0, 1.0) as f32;
        x[i] -= meds[b0] * (1.0 - f) + meds[b1] * f;
    }
}

pub fn robust_std(x: &[f32]) -> f32 {
    let step = (x.len() / 65536).max(1);
    let mut v: Vec<f32> = x.iter().step_by(step).map(|v| v.abs()).collect();
    let s = 1.4826 * median(&mut v);
    if s > 0.0 { s } else { 1.0 }
}

pub fn widths(max_width: usize, snr_loss: f64) -> Vec<usize> {
    let mut out = vec![1usize];
    let f = 1.0 / (1.0 - snr_loss.clamp(0.01, 0.9)).powi(2);
    while *out.last().unwrap() < max_width {
        let last = *out.last().unwrap();
        let next = ((last as f64 * f).ceil() as usize).max(last + 1).min(max_width);
        out.push(next);
    }
    out
}

#[derive(Debug, Clone, Copy)]
pub struct Event {
    pub dm_index: usize,
    pub sample: usize,
    pub width: usize,
    pub snr: f32,
}

pub fn boxcar_noise(cs: &[f64], w: usize) -> (f64, f64) {
    let n = cs.len() - 1;
    let m = n + 1 - w;
    let step = (m / 16384).max(1);
    let mut y: Vec<f64> = (0..m).step_by(step).map(|i| cs[i + w] - cs[i]).collect();
    if y.len() < 32 {
        return (0.0, 1.0 / (w as f64).sqrt());
    }
    let mid = y.len() / 2;
    y.select_nth_unstable_by(mid, |a, b| a.partial_cmp(b).unwrap());
    let med = y[mid];
    let mut d: Vec<f64> = y.iter().map(|v| (v - med).abs()).collect();
    d.select_nth_unstable_by(mid, |a, b| a.partial_cmp(b).unwrap());
    let sd = 1.4826 * d[mid];
    (med, if sd > 0.0 { 1.0 / sd } else { 1.0 / (w as f64).sqrt() })
}

pub fn boxcar_search(x: &[f32], widths: &[usize], threshold: f32) -> Vec<Event> {
    let n = x.len();
    let mut cs = vec![0f64; n + 1];
    for i in 0..n {
        cs[i + 1] = cs[i] + x[i] as f64;
    }
    let mut out = Vec::new();
    for &w in widths {
        if w > n {
            break;
        }
        let (offset, norm) = boxcar_noise(&cs, w);
        let mut best: Option<Event> = None;
        for i in 0..=n - w {
            let s = ((cs[i + w] - cs[i] - offset) * norm) as f32;
            if s >= threshold {
                match best {
                    Some(ref mut b) if i <= b.sample + w => {
                        if s > b.snr {
                            b.sample = i;
                            b.snr = s;
                        }
                    }
                    _ => {
                        if let Some(b) = best.take() {
                            out.push(b);
                        }
                        best = Some(Event { dm_index: 0, sample: i, width: w, snr: s });
                    }
                }
            }
        }
        if let Some(b) = best {
            out.push(b);
        }
    }
    out
}

#[derive(Debug, Clone, Copy)]
pub struct SearchParams {
    pub threshold: f32,
    pub max_width: usize,
    pub snr_loss: f64,
    pub baseline_window: usize,
    pub time_tol: usize,
    pub dm_tol: usize,
    pub width_tol: usize,
}

impl Default for SearchParams {
    fn default() -> Self {
        Self { threshold: 6.0, max_width: 4096, snr_loss: 0.1, baseline_window: 4096, time_tol: 3, dm_tol: 20, width_tol: 3 }
    }
}

pub fn search_trials(trials: &[Vec<f32>], params: &SearchParams) -> Vec<Event> {
    let ws = widths(params.max_width, params.snr_loss);
    trials
        .par_iter()
        .enumerate()
        .flat_map(|(d, t)| {
            let mut x = t.clone();
            baseline_subtract(&mut x, params.baseline_window);
            let s = robust_std(&x);
            for v in x.iter_mut() {
                *v /= s;
            }
            let mut ev = boxcar_search(&x, &ws, params.threshold);
            for e in ev.iter_mut() {
                e.dm_index = d;
            }
            ev
        })
        .collect()
}

pub fn cluster(events: &[Event], params: &SearchParams, widths_list: &[usize]) -> Vec<(Event, u32)> {
    let mut ev = events.to_vec();
    ev.sort_by_key(|e| e.sample);
    let n = ev.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], i: usize) -> usize {
        let mut r = i;
        while p[r] != r {
            p[r] = p[p[r]];
            r = p[r];
        }
        r
    }
    let widx = |w: usize| widths_list.iter().position(|x| *x == w).unwrap_or(0) as isize;
    for i in 0..n {
        let mut j = i + 1;
        while j < n {
            let wmax = ev[i].width.max(ev[j].width);
            if ev[j].sample > ev[i].sample + params.time_tol * wmax {
                break;
            }
            if (ev[i].dm_index as isize - ev[j].dm_index as isize).unsigned_abs() <= params.dm_tol
                && (widx(ev[i].width) - widx(ev[j].width)).unsigned_abs() <= params.width_tol
            {
                let a = find(&mut parent, i);
                let b = find(&mut parent, j);
                if a != b {
                    parent[a] = b;
                }
            }
            j += 1;
        }
    }
    let mut best: std::collections::HashMap<usize, (Event, u32)> = std::collections::HashMap::new();
    for i in 0..n {
        let r = find(&mut parent, i);
        let e = best.entry(r).or_insert((ev[i], 0));
        e.1 += 1;
        if ev[i].snr > e.0.snr {
            e.0 = ev[i];
        }
    }
    let mut out: Vec<(Event, u32)> = best.into_values().collect();
    out.sort_by_key(|(e, _)| e.sample);
    out
}

pub fn to_candidates(clusters: &[(Event, u32)], header: &ObsHeader, start: &Epoch, start_sample: u64, dms: &[f64], f_ref: f64) -> Vec<SinglePulseCandidate> {
    clusters
        .iter()
        .map(|(e, m)| {
            let centre = e.sample as f64 + e.width as f64 / 2.0;
            SinglePulseCandidate {
                obs_id: header.observation_id(),
                dm: dms[e.dm_index],
                epoch: start.add_seconds(centre * header.tsamp),
                f_ref,
                sample: start_sample + centre.floor() as u64,
                width_samples: e.width as u32,
                width_s: e.width as f64 * header.tsamp,
                snr: e.snr as f64,
                dm_index: e.dm_index as u32,
                members: *m,
                beam: header.beam,
                label: None,
                probability: None,
            }
        })
        .collect()
}
