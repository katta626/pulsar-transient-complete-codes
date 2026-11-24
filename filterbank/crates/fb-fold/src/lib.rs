use fb_defs::{Block, Epoch, FbError, Result, Series, K_DM_PULSAR};
use rayon::prelude::*;

#[derive(Debug, Clone, Copy)]
pub struct FoldParams {
    pub f: f64,
    pub fdot: f64,
    pub fddot: f64,
    pub epoch: Epoch,
    pub dm: f64,
    pub nbins: usize,
}

impl FoldParams {
    pub fn phase(&self, t: &Epoch) -> f64 {
        let dt = t.seconds_since(&self.epoch);
        self.f * dt + 0.5 * self.fdot * dt * dt + self.fddot * dt * dt * dt / 6.0
    }
}

#[derive(Debug, Clone)]
pub struct Folded {
    pub nsub: usize,
    pub nchan: usize,
    pub nbins: usize,
    pub data: Vec<f64>,
    pub counts: Vec<u32>,
}

impl Folded {
    pub fn profile(&self) -> Vec<f64> {
        let mut p = vec![0f64; self.nbins];
        for s in 0..self.nsub {
            for c in 0..self.nchan {
                for b in 0..self.nbins {
                    p[b] += self.data[(s * self.nchan + c) * self.nbins + b];
                }
            }
        }
        p
    }

    pub fn reduced_chi2(&self, var_per_sample: f64) -> f64 {
        let p = self.profile();
        let mut cnt = vec![0f64; self.nbins];
        for s in 0..self.nsub {
            for c in 0..self.nchan {
                for b in 0..self.nbins {
                    cnt[b] += self.counts[(s * self.nchan + c) * self.nbins + b] as f64;
                }
            }
        }
        let tot: f64 = p.iter().sum();
        let ntot: f64 = cnt.iter().sum();
        let mean = tot / ntot.max(1.0);
        let chi: f64 = p.iter().zip(&cnt).filter(|(_, n)| **n > 0.0).map(|(v, n)| (v - mean * n).powi(2) / (var_per_sample * n)).sum();
        chi / (self.nbins as f64 - 1.0)
    }
}

pub fn fold_series(s: &Series, fp: &FoldParams, nsub: usize) -> Folded {
    let n = s.len();
    let nsub = nsub.max(1);
    let per = n.div_ceil(nsub);
    let mut data = vec![0f64; nsub * fp.nbins];
    let mut counts = vec![0u32; nsub * fp.nbins];
    for i in 0..n {
        let t = s.header.start.add_seconds(i as f64 * s.tsamp());
        let ph = fp.phase(&t).rem_euclid(1.0);
        let b = ((ph * fp.nbins as f64) as usize).min(fp.nbins - 1);
        let sub = i / per;
        data[sub * fp.nbins + b] += s.data[i] as f64;
        counts[sub * fp.nbins + b] += 1;
    }
    Folded { nsub, nchan: 1, nbins: fp.nbins, data, counts }
}

pub fn fold_block(block: &Block, fp: &FoldParams, nsub: usize, nchan_out: usize, acc: Option<Folded>) -> Result<Folded> {
    if block.npol() != 1 {
        return Err(FbError::Unsupported("fold expects total intensity".into()));
    }
    let h = &block.header;
    let nchan = h.nchan();
    let nchan_out = nchan_out.clamp(1, nchan);
    let per_chan = nchan.div_ceil(nchan_out);
    let f_ref = h.freqs.mhz.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let delays: Vec<f64> = h.freqs.mhz.iter().map(|f| K_DM_PULSAR * fp.dm * (f.powi(-2) - f_ref.powi(-2))).collect();
    let total = h.nsamples.max(1) as usize;
    let per_sub = total.div_ceil(nsub.max(1));
    let mut out = acc.unwrap_or(Folded {
        nsub: nsub.max(1),
        nchan: nchan_out,
        nbins: fp.nbins,
        data: vec![0f64; nsub.max(1) * nchan_out * fp.nbins],
        counts: vec![0u32; nsub.max(1) * nchan_out * fp.nbins],
    });
    let nb = fp.nbins;
    let rows: Vec<(usize, Vec<(usize, f64)>)> = (0..block.nsamp)
        .into_par_iter()
        .map(|t| {
            let gs = block.start_sample as usize + t;
            let sub = (gs / per_sub).min(out.nsub - 1);
            let base = h.epoch_of_sample(gs as u64);
            let spec = block.spectrum(t, 0);
            let v: Vec<(usize, f64)> = (0..nchan)
                .map(|c| {
                    let ph = fp.phase(&base.add_seconds(-delays[c])).rem_euclid(1.0);
                    let b = ((ph * nb as f64) as usize).min(nb - 1);
                    ((c / per_chan) * nb + b, spec[c] as f64)
                })
                .collect();
            (sub, v)
        })
        .collect();
    for (sub, v) in rows {
        for (idx, val) in v {
            let k = sub * nchan_out * nb + idx;
            out.data[k] += val;
            out.counts[k] += 1;
        }
    }
    Ok(out)
}

pub fn normalised_profile(f: &Folded) -> Vec<f64> {
    let p = f.profile();
    let mut cnt = vec![0f64; f.nbins];
    for (i, c) in f.counts.iter().enumerate() {
        cnt[i % f.nbins] += *c as f64;
    }
    p.iter().zip(&cnt).map(|(v, n)| if *n > 0.0 { v / n } else { 0.0 }).collect()
}

pub fn optimise(s: &Series, fp: &FoldParams, df_range: f64, ndf: usize, nfdot: usize, fdot_range: f64, var: f64) -> (FoldParams, f64) {
    let grid: Vec<(f64, f64)> = (0..ndf.max(1))
        .flat_map(|i| (0..nfdot.max(1)).map(move |j| (i, j)))
        .map(|(i, j)| {
            let df = if ndf > 1 { -df_range + 2.0 * df_range * i as f64 / (ndf - 1) as f64 } else { 0.0 };
            let dfd = if nfdot > 1 { -fdot_range + 2.0 * fdot_range * j as f64 / (nfdot - 1) as f64 } else { 0.0 };
            (df, dfd)
        })
        .collect();
    grid.par_iter()
        .map(|(df, dfd)| {
            let p = FoldParams { f: fp.f + df, fdot: fp.fdot + dfd, ..*fp };
            let chi = fold_series(s, &p, 1).reduced_chi2(var);
            (p, chi)
        })
        .reduce(|| (*fp, f64::NEG_INFINITY), |a, b| if b.1 > a.1 { b } else { a })
}
