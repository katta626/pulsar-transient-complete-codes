use fb_defs::{Block, DmPlan, Result};
use rayon::prelude::*;

use crate::{check, DedispOutput};

fn dt_for(max_dt: usize, f_lo: f64, f_hi: f64, f_min: f64, f_max: f64) -> usize {
    (((max_dt as f64 - 1.0) * (f_lo.powi(-2) - f_hi.powi(-2)) / (f_min.powi(-2) - f_max.powi(-2))).ceil().max(0.0)) as usize
}

pub struct Row {
    pub lo: usize,
    pub hi: usize,
    pub shift: usize,
}

pub struct Level {
    pub rows: Vec<Row>,
}

pub struct FdmtPlan {
    pub nf: usize,
    pub dt0: usize,
    pub levels: Vec<Level>,
    pub max_rows: usize,
}

impl FdmtPlan {
    pub fn new(nf: usize, f_min: f64, f_max: f64, max_dt: usize) -> Self {
        let df = (f_max - f_min) / nf as f64;
        let dt0 = dt_for(max_dt, f_min, f_min + df, f_min, f_max);
        let mut ndts: Vec<usize> = vec![dt0 + 1; nf];
        let mut max_rows = nf * (dt0 + 1);
        let mut levels = Vec::new();
        let nlev = nf.trailing_zeros() as usize;
        for it in 1..=nlev {
            let nb = ndts.len() / 2;
            let width = (f_max - f_min) / nb as f64;
            let correction = if it > 1 { df / 2.0 } else { 0.0 };
            let mut prev_off = vec![0usize; ndts.len()];
            for i in 1..ndts.len() {
                prev_off[i] = prev_off[i - 1] + ndts[i - 1];
            }
            let mut rows = Vec::new();
            let mut new_ndts = Vec::with_capacity(nb);
            for b in 0..nb {
                let f_start = f_min + width * b as f64;
                let f_end = f_start + width;
                let f_mid = width / 2.0 + f_start - correction;
                let f_mid_l = width / 2.0 + f_start + correction;
                let ndt = dt_for(max_dt, f_start, f_end, f_min, f_max) + 1;
                let ndt = if it == nlev { ndt.min(max_dt.max(1)) } else { ndt };
                let span = f_end.powi(-2) - f_start.powi(-2);
                for idt in 0..ndt {
                    let dt_mid = ((idt as f64) * (f_mid.powi(-2) - f_start.powi(-2)) / span).round() as usize;
                    let dt_mid_l = ((idt as f64) * (f_mid_l.powi(-2) - f_start.powi(-2)) / span).round() as usize;
                    let dt_rest = idt.saturating_sub(dt_mid_l);
                    rows.push(Row {
                        lo: prev_off[2 * b] + dt_mid.min(ndts[2 * b] - 1),
                        hi: prev_off[2 * b + 1] + dt_rest.min(ndts[2 * b + 1] - 1),
                        shift: dt_mid_l,
                    });
                }
                new_ndts.push(ndt);
            }
            max_rows = max_rows.max(rows.len());
            levels.push(Level { rows });
            ndts = new_ndts;
        }
        Self { nf, dt0, levels, max_rows }
    }

    pub fn ndt(&self) -> usize {
        self.levels.last().map(|l| l.rows.len()).unwrap_or(self.dt0 + 1)
    }

    pub fn execute(&self, image: &[f32], nt: usize) -> Vec<f32> {
        let mut a = vec![0f32; self.max_rows * nt];
        let mut b = vec![0f32; self.max_rows * nt];
        let ndt0 = self.dt0 + 1;
        a[..self.nf * ndt0 * nt].par_chunks_mut(ndt0 * nt).enumerate().for_each(|(c, band)| {
            let row = &image[c * nt..(c + 1) * nt];
            band[..nt].copy_from_slice(row);
            for k in 1..ndt0 {
                let (prev, cur) = band.split_at_mut(k * nt);
                let prev = &prev[(k - 1) * nt..];
                for t in k..nt {
                    cur[t] = prev[t] + row[t - k];
                }
            }
        });
        for lev in &self.levels {
            let src = &a;
            b[..lev.rows.len() * nt].par_chunks_mut(nt).zip(lev.rows.par_iter()).for_each(|(out, r)| {
                let x = &src[r.lo * nt..(r.lo + 1) * nt];
                let y = &src[r.hi * nt..(r.hi + 1) * nt];
                let split = r.shift.min(nt);
                out[..split].copy_from_slice(&x[..split]);
                for ((o, p), q) in out[split..].iter_mut().zip(&x[split..]).zip(&y[..nt - split]) {
                    *o = p + q;
                }
            });
            std::mem::swap(&mut a, &mut b);
        }
        a.truncate(self.ndt() * nt);
        a
    }
}

pub fn transform(image: &[f32], nf: usize, nt: usize, f_min: f64, f_max: f64, max_dt: usize) -> Vec<f32> {
    FdmtPlan::new(nf, f_min, f_max, max_dt).execute(image, nt)
}

pub struct FdmtSetup {
    pub plan: FdmtPlan,
    pub order: Vec<i64>,
    pub out_rows: Vec<usize>,
    pub out_shift: Vec<usize>,
    pub nout: usize,
}

impl FdmtSetup {
    pub fn new(h: &fb_defs::ObsHeader, plan: &DmPlan, nt: usize, killmask: Option<&[bool]>) -> Self {
        let nchan = h.nchan();
        let bw = h.freqs.chan_bw.abs();
        let mut idx: Vec<usize> = (0..nchan).collect();
        idx.sort_by(|a, b| h.freqs.mhz[*a].partial_cmp(&h.freqs.mhz[*b]).unwrap());
        let f_min = h.freqs.mhz[idx[0]] - bw / 2.0;
        let nf = nchan.next_power_of_two();
        let f_max = f_min + nf as f64 * bw;
        let tsamp = plan.tsamp;
        let dm_max = plan.dms.iter().cloned().fold(0.0, f64::max);
        let max_dt = ((plan.k_dm * dm_max * (f_min.powi(-2) - f_max.powi(-2)) / tsamp).ceil() as usize + 2).min(nt);
        let fp = FdmtPlan::new(nf, f_min, f_max, max_dt);
        let ndt = fp.ndt();
        let order = (0..nf)
            .map(|r| if r < nchan && !killmask.map(|k| k[idx[r]]).unwrap_or(false) { idx[r] as i64 } else { -1 })
            .collect();
        let out_rows = plan.dms.iter().map(|dm| ((plan.k_dm * dm * (f_min.powi(-2) - f_max.powi(-2)) / tsamp).round() as usize).min(ndt - 1)).collect();
        let out_shift = plan.dms.iter().map(|dm| (plan.k_dm * dm * (f_min.powi(-2) - plan.f_ref.powi(-2)) / tsamp).round() as usize).collect();
        Self { plan: fp, order, out_rows, out_shift, nout: nt - plan.max_delay_samples() }
    }
}

pub fn fdmt(block: &Block, plan: &DmPlan, killmask: Option<&[bool]>) -> Result<DedispOutput> {
    fdmt_chunked(block, plan, killmask, 1 << 28)
}

pub fn fdmt_chunked(block: &Block, plan: &DmPlan, killmask: Option<&[bool]>, budget_floats: usize) -> Result<DedispOutput> {
    check(block, plan)?;
    let nchan = block.header.nchan();
    let nt = block.nsamp;
    let st = FdmtSetup::new(&block.header, plan, nt, killmask);
    let nf = st.plan.nf;
    let depth = st.plan.ndt() + 2;
    let per_t = 2 * st.plan.max_rows.max(nf);
    let nc = (budget_floats / per_t).saturating_sub(depth).max(depth).max(256).min(nt);
    let nout = st.nout;
    let mut data = vec![0f32; plan.ndm() * nout];
    let mut c0 = 0usize;
    while c0 < nt {
        let c1 = (c0 + nc).min(nt);
        let s = c0 - depth.min(c0);
        let len = c1 - s;
        let mut image = vec![0f32; nf * len];
        image.par_chunks_mut(len).enumerate().for_each(|(row, dst)| {
            if st.order[row] < 0 {
                return;
            }
            let c = st.order[row] as usize;
            for (i, v) in dst.iter_mut().enumerate() {
                *v = block.data[(s + i) * nchan + c];
            }
        });
        let dmt = st.plan.execute(&image, len);
        data.par_chunks_mut(nout).enumerate().for_each(|(d, out)| {
            let row = &dmt[st.out_rows[d] * len..(st.out_rows[d] + 1) * len];
            let shift = st.out_shift[d];
            for tp in c0.max(shift)..c1 {
                let t = tp - shift;
                if t < nout {
                    out[t] = row[tp - s];
                }
            }
        });
        c0 = c1;
    }
    Ok(DedispOutput { ndm: plan.ndm(), nout, data })
}
