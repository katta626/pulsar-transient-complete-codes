use serde::{Deserialize, Serialize};

use crate::consts::K_DM_PULSAR;
use crate::error::{FbError, Result};
use crate::header::FreqAxis;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DmPlan {
    pub freqs: FreqAxis,
    pub tsamp: f64,
    pub f_ref: f64,
    pub dms: Vec<f64>,
    pub k_dm: f64,
}

impl DmPlan {
    pub fn from_list(freqs: FreqAxis, tsamp: f64, dms: Vec<f64>) -> Result<Self> {
        if dms.is_empty() {
            return Err(FbError::Invalid("empty DM list".into()));
        }
        let f_ref = freqs.mhz.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        Ok(Self { freqs, tsamp, f_ref, dms, k_dm: K_DM_PULSAR })
    }

    pub fn linear(freqs: FreqAxis, tsamp: f64, lo: f64, hi: f64, step: f64) -> Result<Self> {
        if step <= 0.0 || hi < lo {
            return Err(FbError::Invalid("DM range must satisfy lo <= hi and step > 0".into()));
        }
        let n = ((hi - lo) / step).floor() as usize + 1;
        Self::from_list(freqs, tsamp, (0..n).map(|i| lo + step * i as f64).collect())
    }

    pub fn tolerance(freqs: FreqAxis, tsamp: f64, lo: f64, hi: f64, pulse_width: f64, tol: f64) -> Result<Self> {
        if tol <= 1.0 {
            return Err(FbError::Invalid("tolerance must exceed 1 (allowed S/N loss factor)".into()));
        }
        let f_lo = freqs.mhz.iter().cloned().fold(f64::INFINITY, f64::min);
        let f_hi = freqs.mhz.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let f_c = 0.5 * (f_lo + f_hi);
        let chan_bw = freqs.chan_bw.abs();
        let sweep = K_DM_PULSAR * (f_lo.powi(-2) - f_hi.powi(-2));
        let mut dms = vec![lo];
        let mut dm = lo;
        while dm < hi {
            let t_chan = 2.0 * K_DM_PULSAR * dm * chan_bw / f_c.powi(3);
            let w_eff = (tsamp * tsamp + pulse_width * pulse_width + t_chan * t_chan).sqrt();
            let step = 2.0 * w_eff * (tol * tol - 1.0).sqrt() / sweep;
            dm += step.max(1e-6);
            dms.push(dm.min(hi));
        }
        dms.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
        Self::from_list(freqs, tsamp, dms)
    }

    pub fn ndm(&self) -> usize {
        self.dms.len()
    }

    pub fn delay_seconds(&self, dm: f64, chan: usize) -> f64 {
        self.k_dm * dm * (self.freqs.mhz[chan].powi(-2) - self.f_ref.powi(-2))
    }

    pub fn delay_table(&self) -> Vec<i32> {
        let nchan = self.freqs.nchan();
        let mut out = Vec::with_capacity(self.ndm() * nchan);
        for &dm in &self.dms {
            for c in 0..nchan {
                out.push((self.delay_seconds(dm, c) / self.tsamp).round() as i32);
            }
        }
        out
    }

    pub fn max_delay_samples(&self) -> usize {
        let dm = self.dms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        (0..self.freqs.nchan()).map(|c| (self.delay_seconds(dm, c) / self.tsamp).round() as i64).max().unwrap_or(0).max(0) as usize
    }
}
