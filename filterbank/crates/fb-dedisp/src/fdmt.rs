use fb_defs::{Block, DmPlan, Result};
use rayon::prelude::*;

use crate::{check, DedispOutput};

struct Band {
    ndt: usize,
    data: Vec<f32>,
}

fn dt_for(max_dt: usize, f_lo: f64, f_hi: f64, f_min: f64, f_max: f64) -> usize {
    (((max_dt as f64 - 1.0) * (f_lo.powi(-2) - f_hi.powi(-2)) / (f_min.powi(-2) - f_max.powi(-2))).ceil().max(0.0)) as usize
}

pub fn transform(image: &[f32], nf: usize, nt: usize, f_min: f64, f_max: f64, max_dt: usize) -> Vec<f32> {
    let df = (f_max - f_min) / nf as f64;
    let dt0 = dt_for(max_dt, f_min, f_min + df, f_min, f_max);
    let mut bands: Vec<Band> = (0..nf)
        .into_par_iter()
        .map(|c| {
            let row = &image[c * nt..(c + 1) * nt];
            let ndt = dt0 + 1;
            let mut data = vec![0f32; ndt * nt];
            data[..nt].copy_from_slice(row);
            for k in 1..ndt {
                let (prev, cur) = data.split_at_mut(k * nt);
                let prev = &prev[(k - 1) * nt..];
                for t in k..nt {
                    cur[t] = prev[t] + row[t - k];
                }
            }
            Band { ndt, data }
        })
        .collect();
    let levels = nf.trailing_zeros() as usize;
    for it in 1..=levels {
        let nb = bands.len() / 2;
        let width = (f_max - f_min) / nb as f64;
        let correction = if it > 1 { df / 2.0 } else { 0.0 };
        let next: Vec<Band> = (0..nb)
            .into_par_iter()
            .map(|b| {
                let lo = &bands[2 * b];
                let hi = &bands[2 * b + 1];
                let f_start = f_min + width * b as f64;
                let f_end = f_start + width;
                let f_mid = (f_end - f_start) / 2.0 + f_start - correction;
                let f_mid_l = (f_end - f_start) / 2.0 + f_start + correction;
                let local = dt_for(max_dt, f_start, f_end, f_min, f_max);
                let ndt = local + 1;
                let mut data = vec![0f32; ndt * nt];
                let span = f_end.powi(-2) - f_start.powi(-2);
                for idt in 0..ndt {
                    let dt_mid = ((idt as f64) * (f_mid.powi(-2) - f_start.powi(-2)) / span).round() as usize;
                    let dt_mid_l = ((idt as f64) * (f_mid_l.powi(-2) - f_start.powi(-2)) / span).round() as usize;
                    let dt_rest = idt.saturating_sub(dt_mid_l);
                    let a = &lo.data[dt_mid.min(lo.ndt - 1) * nt..(dt_mid.min(lo.ndt - 1) + 1) * nt];
                    let bh = &hi.data[dt_rest.min(hi.ndt - 1) * nt..(dt_rest.min(hi.ndt - 1) + 1) * nt];
                    let out = &mut data[idt * nt..(idt + 1) * nt];
                    let split = dt_mid_l.min(nt);
                    out[..split].copy_from_slice(&a[..split]);
                    for t in split..nt {
                        out[t] = a[t] + bh[t - dt_mid_l];
                    }
                }
                Band { ndt, data }
            })
            .collect();
        bands = next;
    }
    let b = bands.pop().unwrap();
    let mut out = b.data;
    out.truncate(max_dt.min(b.ndt) * nt);
    out
}

pub fn fdmt(block: &Block, plan: &DmPlan, killmask: Option<&[bool]>) -> Result<DedispOutput> {
    check(block, plan)?;
    let h = &block.header;
    let nchan = h.nchan();
    let nt = block.nsamp;
    let bw = h.freqs.chan_bw.abs();
    let ascending: Vec<usize> = {
        let mut idx: Vec<usize> = (0..nchan).collect();
        idx.sort_by(|a, b| h.freqs.mhz[*a].partial_cmp(&h.freqs.mhz[*b]).unwrap());
        idx
    };
    let f_min = h.freqs.mhz[ascending[0]] - bw / 2.0;
    let nf = nchan.next_power_of_two();
    let f_max = f_min + nf as f64 * bw;
    let tsamp = plan.tsamp;
    let dm_max = plan.dms.iter().cloned().fold(0.0, f64::max);
    let max_dt = ((plan.k_dm * dm_max * (f_min.powi(-2) - f_max.powi(-2)) / tsamp).ceil() as usize + 2).min(nt);
    let mut image = vec![0f32; nf * nt];
    image.par_chunks_mut(nt).enumerate().for_each(|(row, dst)| {
        if row >= nchan {
            return;
        }
        let c = ascending[row];
        if killmask.map(|k| k[c]).unwrap_or(false) {
            return;
        }
        for (t, v) in dst.iter_mut().enumerate() {
            *v = block.data[t * nchan + c];
        }
    });
    let dmt = transform(&image, nf, nt, f_min, f_max, max_dt);
    let ndt = dmt.len() / nt;
    let maxd = plan.max_delay_samples();
    let nout = nt - maxd;
    let mut data = vec![0f32; plan.ndm() * nout];
    data.par_chunks_mut(nout).enumerate().for_each(|(d, out)| {
        let dm = plan.dms[d];
        let dt = ((plan.k_dm * dm * (f_min.powi(-2) - f_max.powi(-2)) / tsamp).round() as usize).min(ndt - 1);
        let shift = (plan.k_dm * dm * (f_min.powi(-2) - plan.f_ref.powi(-2)) / tsamp).round() as usize;
        let row = &dmt[dt * nt..(dt + 1) * nt];
        for (t, o) in out.iter_mut().enumerate() {
            let src = t + shift;
            *o = if src < nt { row[src] } else { 0.0 };
        }
    });
    Ok(DedispOutput { ndm: plan.ndm(), nout, data })
}
