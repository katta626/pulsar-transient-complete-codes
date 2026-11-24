pub mod fdmt;

pub use fdmt::fdmt;

use fb_defs::{Block, DmPlan, FbError, Result, Series};
use rayon::prelude::*;

pub struct DedispOutput {
    pub ndm: usize,
    pub nout: usize,
    pub data: Vec<f32>,
}

impl DedispOutput {
    pub fn trial(&self, i: usize) -> &[f32] {
        &self.data[i * self.nout..(i + 1) * self.nout]
    }
}

pub fn check(block: &Block, plan: &DmPlan) -> Result<()> {
    if block.npol() != 1 {
        return Err(FbError::Unsupported("dedispersion expects total intensity (npol=1)".into()));
    }
    block.header.check_same_freqs(&plan.freqs)?;
    if (block.header.tsamp - plan.tsamp).abs() > 1e-12 * plan.tsamp.max(1.0) {
        return Err(FbError::Incompatible(format!("plan tsamp {} != data tsamp {}", plan.tsamp, block.header.tsamp)));
    }
    if block.nsamp <= plan.max_delay_samples() {
        return Err(FbError::Invalid(format!(
            "block of {} samples is shorter than the maximum delay {}",
            block.nsamp,
            plan.max_delay_samples()
        )));
    }
    Ok(())
}

pub fn transpose(block: &Block) -> Vec<f32> {
    let nchan = block.nchan();
    let nsamp = block.nsamp;
    let mut out = vec![0f32; nchan * nsamp];
    out.par_chunks_mut(nsamp).enumerate().for_each(|(c, col)| {
        for (t, v) in col.iter_mut().enumerate() {
            *v = block.data[t * nchan + c];
        }
    });
    out
}

pub fn brute_force(block: &Block, plan: &DmPlan, killmask: Option<&[bool]>) -> Result<DedispOutput> {
    check(block, plan)?;
    let nchan = block.nchan();
    let nsamp = block.nsamp;
    let maxd = plan.max_delay_samples();
    let nout = nsamp - maxd;
    let delays = plan.delay_table();
    let tr = transpose(block);
    let use_chan: Vec<bool> = match killmask {
        Some(k) => k.iter().map(|b| !b).collect(),
        None => vec![true; nchan],
    };
    let mut data = vec![0f32; plan.ndm() * nout];
    data.par_chunks_mut(nout).enumerate().for_each(|(d, out)| {
        let row = &delays[d * nchan..(d + 1) * nchan];
        for c in 0..nchan {
            if !use_chan[c] {
                continue;
            }
            let off = row[c].max(0) as usize;
            let src = &tr[c * nsamp + off..c * nsamp + off + nout];
            for (o, s) in out.iter_mut().zip(src) {
                *o += *s;
            }
        }
    });
    Ok(DedispOutput { ndm: plan.ndm(), nout, data })
}

pub fn subband(block: &Block, plan: &DmPlan, nsub: usize, killmask: Option<&[bool]>) -> Result<DedispOutput> {
    check(block, plan)?;
    let nchan = block.nchan();
    let nsamp = block.nsamp;
    let nsub = nsub.clamp(1, nchan);
    let per = nchan.div_ceil(nsub);
    let maxd = plan.max_delay_samples();
    let nout = nsamp - maxd;
    let tsamp = plan.tsamp;
    let tr = transpose(block);
    let use_chan: Vec<bool> = match killmask {
        Some(k) => k.iter().map(|b| !b).collect(),
        None => vec![true; nchan],
    };
    let groups: Vec<(usize, usize)> = (0..nsub).map(|s| (s * per, ((s + 1) * per).min(nchan))).filter(|g| g.0 < g.1).collect();
    let ref_chan: Vec<usize> = groups
        .iter()
        .map(|g| (g.0..g.1).max_by(|a, b| plan.freqs.mhz[*a].partial_cmp(&plan.freqs.mhz[*b]).unwrap()).unwrap())
        .collect();
    let ndm = plan.ndm();
    let step = ((ndm as f64).sqrt().ceil() as usize).max(1);
    let mut data = vec![0f32; ndm * nout];
    let blocks: Vec<(usize, usize)> = (0..ndm).step_by(step).map(|i| (i, (i + step).min(ndm))).collect();
    for (d0, d1) in blocks {
        let g_dm = plan.dms[(d0 + d1) / 2];
        let sub_len = nsamp;
        let subs: Vec<Vec<f32>> = groups
            .par_iter()
            .enumerate()
            .map(|(gi, g)| {
                let rc = ref_chan[gi];
                let mut acc = vec![0f32; sub_len];
                for c in g.0..g.1 {
                    if !use_chan[c] {
                        continue;
                    }
                    let fine = ((plan.delay_seconds(g_dm, c) - plan.delay_seconds(g_dm, rc)) / tsamp).round().max(0.0) as usize;
                    let src = &tr[c * nsamp..(c + 1) * nsamp];
                    let n = nsamp - fine;
                    for (o, s) in acc[..n].iter_mut().zip(&src[fine..]) {
                        *o += *s;
                    }
                }
                acc
            })
            .collect();
        data[d0 * nout..d1 * nout].par_chunks_mut(nout).enumerate().for_each(|(k, out)| {
            let d = d0 + k;
            for (gi, sub) in subs.iter().enumerate() {
                let coarse = ((plan.delay_seconds(plan.dms[d], ref_chan[gi])) / tsamp).round().clamp(0.0, maxd as f64) as usize;
                for (o, s) in out.iter_mut().zip(&sub[coarse..coarse + nout]) {
                    *o += *s;
                }
            }
        });
    }
    Ok(DedispOutput { ndm, nout, data })
}

pub fn to_series(block: &Block, plan: &DmPlan, out: &DedispOutput, i: usize) -> Series {
    let mut h = (*block.header).clone();
    h.nsamples = out.nout as u64;
    h.start = block.start_epoch();
    h.provenance.push(format!("dedisp dm={:.4}", plan.dms[i]));
    Series { header: h, dm: plan.dms[i], f_ref: plan.f_ref, data: out.trial(i).to_vec() }
}
