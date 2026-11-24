pub mod iqrm;
pub mod stats;

use fb_defs::{Block, FbError, ObsHeader, Result, RfiMask};
use rayon::prelude::*;

pub use iqrm::iqrm;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Statistic {
    Std,
    Mean,
    Skew,
    Kurtosis,
}

fn intensity(block: &Block) -> Result<()> {
    if block.npol() != 1 {
        return Err(FbError::Unsupported("RFI stages expect total-intensity blocks (npol=1)".into()));
    }
    Ok(())
}

pub fn zero_dm(block: &mut Block, keep_channel_mean: bool) -> Result<()> {
    intensity(block)?;
    let nchan = block.nchan();
    let nsamp = block.nsamp;
    let means: Vec<f32> = if keep_channel_mean {
        let st = stats::channel_stats(&block.data, nsamp, nchan);
        st.mean.iter().map(|m| *m as f32).collect()
    } else {
        vec![0.0; nchan]
    };
    block.data.par_chunks_mut(nchan).for_each(|row| {
        let m = row.iter().sum::<f32>() / nchan as f32;
        for (v, cm) in row.iter_mut().zip(&means) {
            *v = *v - m + *cm;
        }
    });
    let mut h = (*block.header).clone();
    h.provenance.push("rfi:zero_dm".into());
    block.header = std::sync::Arc::new(h);
    Ok(())
}

pub fn spectral_kurtosis(block: &Block, m: usize, n_integrated: f64, sigma: f64, mask: &mut RfiMask) -> Result<()> {
    intensity(block)?;
    mask.check_applies_to(&block.header)?;
    let nchan = block.nchan();
    let nwin = block.nsamp / m.max(2);
    let mf = m as f64;
    let nd = n_integrated.max(1.0);
    let var = 2.0 * nd * (nd + 1.0) * mf * mf / ((mf - 1.0) * (6.0 + 5.0 * mf * nd + mf * mf * nd * nd));
    let bound = sigma * var.sqrt();
    for w in 0..nwin {
        let t0 = w * m;
        let mut s1 = vec![0f64; nchan];
        let mut s2 = vec![0f64; nchan];
        for t in t0..t0 + m {
            let row = block.spectrum(t, 0);
            for c in 0..nchan {
                let x = row[c] as f64;
                s1[c] += x;
                s2[c] += x * x;
            }
        }
        let t_mid = block.start_sample + (t0 + m / 2) as u64;
        let epoch = block.header.epoch_of_sample(t_mid);
        let int = (epoch.seconds_since(&mask.start) / mask.interval).floor();
        if int < 0.0 || int as usize >= mask.nint {
            continue;
        }
        for c in 0..nchan {
            if s1[c] <= 0.0 {
                continue;
            }
            let sk = (mf * nd + 1.0) / (mf - 1.0) * (mf * s2[c] / (s1[c] * s1[c]) - 1.0);
            if (sk - 1.0).abs() > bound {
                let (lo, hi) = block.header.freqs.edges(c);
                for mc in 0..mask.nchan() {
                    if mask.chan_hi[mc] > lo && mask.chan_lo[mc] < hi {
                        mask.set(int as usize, mc, true);
                    }
                }
            }
        }
    }
    mask.provenance.push(format!("rfi:sk m={m} n={nd} sigma={sigma}"));
    Ok(())
}

pub fn iqrm_mask(block: &Block, stat: Statistic, radius: usize, threshold: f64, mask: &mut RfiMask) -> Result<usize> {
    intensity(block)?;
    mask.check_applies_to(&block.header)?;
    let nchan = block.nchan();
    let st = stats::channel_stats(&block.data, block.nsamp, nchan);
    let x = match stat {
        Statistic::Std => st.std,
        Statistic::Mean => st.mean,
        Statistic::Skew => st.skew.iter().map(|v| v.abs()).collect(),
        Statistic::Kurtosis => st.kurt.iter().map(|v| v.abs()).collect(),
    };
    let bad = iqrm(&x, radius, threshold);
    let t0 = block.start_epoch();
    let t1 = block.header.epoch_of_sample(block.start_sample + block.nsamp as u64);
    let a = (t0.seconds_since(&mask.start) / mask.interval).floor().max(0.0) as usize;
    let b = ((t1.seconds_since(&mask.start) / mask.interval).ceil().max(0.0) as usize).min(mask.nint);
    let mut n = 0;
    for (c, flag) in bad.iter().enumerate() {
        if !flag {
            continue;
        }
        n += 1;
        let (lo, hi) = block.header.freqs.edges(c);
        for mc in 0..mask.nchan() {
            if mask.chan_hi[mc] > lo && mask.chan_lo[mc] < hi {
                for i in a..b {
                    mask.set(i, mc, true);
                }
            }
        }
    }
    mask.provenance.push(format!("rfi:iqrm stat={stat:?} radius={radius} t={threshold}"));
    Ok(n)
}

pub fn interval_stats_mask(block: &Block, threshold: f64, mask: &mut RfiMask) -> Result<()> {
    intensity(block)?;
    mask.check_applies_to(&block.header)?;
    let nchan = block.nchan();
    let per = ((mask.interval / block.header.tsamp).round() as usize).max(2);
    let nwin = block.nsamp.div_ceil(per);
    let mut means = vec![0f32; nwin * nchan];
    let mut stds = vec![0f32; nwin * nchan];
    for w in 0..nwin {
        let t0 = w * per;
        let t1 = (t0 + per).min(block.nsamp);
        let n = (t1 - t0) as f64;
        for c in 0..nchan {
            let mut s1 = 0f64;
            let mut s2 = 0f64;
            for t in t0..t1 {
                let x = block.data[t * nchan + c] as f64;
                s1 += x;
                s2 += x * x;
            }
            let m = s1 / n;
            means[w * nchan + c] = m as f32;
            stds[w * nchan + c] = ((s2 / n - m * m).max(0.0)).sqrt() as f32;
        }
    }
    let (mm, ms) = stats::median_mad(&means);
    let (sm, ss) = stats::median_mad(&stds);
    for w in 0..nwin {
        let t_mid = block.start_sample + (w * per + per / 2).min(block.nsamp.saturating_sub(1)) as u64;
        let epoch = block.header.epoch_of_sample(t_mid);
        let int = (epoch.seconds_since(&mask.start) / mask.interval).floor();
        if int < 0.0 || int as usize >= mask.nint {
            continue;
        }
        for c in 0..nchan {
            let bm = ms > 0.0 && ((means[w * nchan + c] - mm).abs() / ms) as f64 > threshold;
            let bs = ss > 0.0 && ((stds[w * nchan + c] - sm).abs() / ss) as f64 > threshold;
            if bm || bs {
                let (lo, hi) = block.header.freqs.edges(c);
                for mc in 0..mask.nchan() {
                    if mask.chan_hi[mc] > lo && mask.chan_lo[mc] < hi {
                        mask.set(int as usize, mc, true);
                    }
                }
            }
        }
    }
    mask.provenance.push(format!("rfi:interval_stats t={threshold}"));
    Ok(())
}

pub fn normalise(block: &mut Block) -> Result<(Vec<f32>, Vec<f32>)> {
    intensity(block)?;
    let nchan = block.nchan();
    let nsamp = block.nsamp;
    let cols: Vec<(f32, f32)> = (0..nchan)
        .into_par_iter()
        .map(|c| {
            let v: Vec<f32> = (0..nsamp).map(|t| block.data[t * nchan + c]).collect();
            let (m, s) = stats::median_mad(&v);
            (m, if s > 0.0 { s } else { 1.0 })
        })
        .collect();
    let (med, sd): (Vec<f32>, Vec<f32>) = cols.into_iter().unzip();
    block.data.par_chunks_mut(nchan).for_each(|row| {
        for c in 0..nchan {
            row[c] = (row[c] - med[c]) / sd[c];
        }
    });
    let mut h = (*block.header).clone();
    h.provenance.push("rfi:normalise median/mad".into());
    block.header = std::sync::Arc::new(h);
    Ok((med, sd))
}

pub fn new_mask(header: &ObsHeader, interval: f64) -> RfiMask {
    RfiMask::empty(header, interval, "filterbank")
}
