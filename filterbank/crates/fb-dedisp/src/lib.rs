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
    let tr = transpose(block);
    let use_chan: Vec<bool> = match killmask {
        Some(k) => k.iter().map(|b| !b).collect(),
        None => vec![true; nchan],
    };
    let mut data = vec![0f32; plan.ndm() * nout];
    let groups: Vec<(usize, usize)> = (0..nsub).map(|s| (s * per, ((s + 1) * per).min(nchan))).filter(|g| g.0 < g.1).collect();
    let ref_chan: Vec<usize> = groups.iter().map(|g| (g.0 + g.1) / 2).collect();
    let nominal: Vec<f64> = {
        let mut v = Vec::new();
        let step = ((plan.ndm() as f64).sqrt().ceil() as usize).max(1);
        let mut i = 0;
        while i < plan.ndm() {
            v.push(plan.dms[(i + step / 2).min(plan.ndm() - 1)]);
            i += step;
        }
        v
    };
    let step = ((plan.ndm() as f64).sqrt().ceil() as usize).max(1);
    let tsamp = plan.tsamp;
    data.par_chunks_mut(nout).enumerate().for_each(|(d, out)| {
        let g_dm = nominal[(d / step).min(nominal.len() - 1)];
        for (gi, g) in groups.iter().enumerate() {
            let rc = ref_chan[gi];
            let coarse = (plan.delay_seconds(plan.dms[d], rc) / tsamp).round() as i64;
            for c in g.0..g.1 {
                if !use_chan[c] {
                    continue;
                }
                let fine = ((plan.delay_seconds(g_dm, c) - plan.delay_seconds(g_dm, rc)) / tsamp).round() as i64;
                let off = (coarse + fine).clamp(0, maxd as i64) as usize;
                let src = &tr[c * nsamp + off..c * nsamp + off + nout];
                for (o, s) in out.iter_mut().zip(src) {
                    *o += *s;
                }
            }
        }
    });
    Ok(DedispOutput { ndm: plan.ndm(), nout, data })
}

pub fn to_series(block: &Block, plan: &DmPlan, out: &DedispOutput, i: usize) -> Series {
    let mut h = (*block.header).clone();
    h.nsamples = out.nout as u64;
    h.start = block.start_epoch();
    h.provenance.push(format!("dedisp dm={:.4}", plan.dms[i]));
    Series { header: h, dm: plan.dms[i], f_ref: plan.f_ref, data: out.trial(i).to_vec() }
}
