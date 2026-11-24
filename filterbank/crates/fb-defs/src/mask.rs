use serde::{Deserialize, Serialize};

use crate::error::{FbError, Result};
use crate::header::{FreqAxis, ObsHeader};
use crate::time::Epoch;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RfiMask {
    pub obs_id: String,
    pub start: Epoch,
    pub interval: f64,
    pub nint: usize,
    pub chan_lo: Vec<f64>,
    pub chan_hi: Vec<f64>,
    pub bad: Vec<bool>,
    pub method: String,
    pub provenance: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct MaskSegment {
    pub t0: usize,
    pub t1: usize,
    pub bad: Vec<bool>,
}

impl RfiMask {
    pub fn empty(header: &ObsHeader, interval: f64, method: &str) -> Self {
        let nint = ((header.duration() / interval).ceil() as usize).max(1);
        let n = header.nchan();
        let (chan_lo, chan_hi): (Vec<f64>, Vec<f64>) = (0..n).map(|c| header.freqs.edges(c)).unzip();
        Self {
            obs_id: header.observation_id(),
            start: header.start,
            interval,
            nint,
            chan_lo,
            chan_hi,
            bad: vec![false; nint * n],
            method: method.into(),
            provenance: header.provenance.clone(),
        }
    }

    pub fn from_freqs(obs_id: String, start: Epoch, interval: f64, nint: usize, freqs: &FreqAxis, method: &str) -> Self {
        let n = freqs.nchan();
        let (chan_lo, chan_hi): (Vec<f64>, Vec<f64>) = (0..n).map(|c| freqs.edges(c)).unzip();
        Self { obs_id, start, interval, nint, chan_lo, chan_hi, bad: vec![false; nint * n], method: method.into(), provenance: vec![] }
    }

    pub fn nchan(&self) -> usize {
        self.chan_lo.len()
    }

    pub fn set(&mut self, int: usize, chan: usize, bad: bool) {
        let n = self.nchan();
        self.bad[int * n + chan] = bad;
    }

    pub fn get(&self, int: usize, chan: usize) -> bool {
        self.bad[int * self.nchan() + chan]
    }

    pub fn zap_channel(&mut self, chan: usize) {
        for i in 0..self.nint {
            self.set(i, chan, true);
        }
    }

    pub fn zap_interval(&mut self, int: usize) {
        for c in 0..self.nchan() {
            self.set(int, c, true);
        }
    }

    pub fn zap_freq_range(&mut self, lo_mhz: f64, hi_mhz: f64) {
        for c in 0..self.nchan() {
            if self.chan_hi[c] > lo_mhz && self.chan_lo[c] < hi_mhz {
                self.zap_channel(c);
            }
        }
    }

    pub fn fraction_bad(&self) -> f64 {
        self.bad.iter().filter(|b| **b).count() as f64 / self.bad.len().max(1) as f64
    }

    pub fn union(&self, other: &RfiMask) -> Result<RfiMask> {
        if self.obs_id != other.obs_id {
            return Err(FbError::Incompatible(format!("masks belong to '{}' and '{}'", self.obs_id, other.obs_id)));
        }
        if self.nchan() == other.nchan() && self.nint == other.nint && (self.interval - other.interval).abs() < 1e-9 {
            let mut m = self.clone();
            for (a, b) in m.bad.iter_mut().zip(&other.bad) {
                *a |= *b;
            }
            m.method = format!("{}+{}", self.method, other.method);
            return Ok(m);
        }
        let mut m = self.clone();
        for i in 0..self.nint {
            let t_lo = self.start.add_seconds(i as f64 * self.interval);
            let t_hi = t_lo.add_seconds(self.interval);
            for c in 0..self.nchan() {
                if other.is_bad_region(&t_lo, &t_hi, self.chan_lo[c], self.chan_hi[c]) {
                    m.set(i, c, true);
                }
            }
        }
        m.method = format!("{}+{}", self.method, other.method);
        Ok(m)
    }

    pub fn is_bad_region(&self, t_lo: &Epoch, t_hi: &Epoch, f_lo: f64, f_hi: f64) -> bool {
        let a = (t_lo.seconds_since(&self.start) / self.interval).floor().max(0.0) as usize;
        let b = ((t_hi.seconds_since(&self.start) / self.interval).ceil().max(0.0) as usize).min(self.nint);
        let chans: Vec<usize> = (0..self.nchan()).filter(|&c| self.chan_hi[c] > f_lo && self.chan_lo[c] < f_hi).collect();
        (a..b).any(|i| chans.iter().any(|&c| self.get(i, c)))
    }

    pub fn check_applies_to(&self, header: &ObsHeader) -> Result<()> {
        if self.obs_id != header.observation_id() {
            return Err(FbError::Incompatible(format!(
                "mask was made for '{}' but data are '{}'",
                self.obs_id,
                header.observation_id()
            )));
        }
        let m_lo = self.chan_lo.iter().cloned().fold(f64::INFINITY, f64::min);
        let m_hi = self.chan_hi.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        if header.freqs.hi() <= m_lo || header.freqs.lo() >= m_hi {
            return Err(FbError::Incompatible("mask and data frequency ranges do not overlap".into()));
        }
        Ok(())
    }

    pub fn to_indices(&self, header: &ObsHeader, start_sample: u64, nsamp: usize) -> Result<Vec<MaskSegment>> {
        let nchan = header.nchan();
        let map: Vec<Vec<usize>> = (0..nchan)
            .map(|c| {
                let (lo, hi) = header.freqs.edges(c);
                (0..self.nchan()).filter(|&m| self.chan_hi[m] > lo && self.chan_lo[m] < hi).collect()
            })
            .collect();
        let mut out: Vec<MaskSegment> = Vec::new();
        let mut t = 0usize;
        while t < nsamp {
            let epoch = header.epoch_of_sample(start_sample + t as u64);
            let off = epoch.seconds_since(&self.start);
            let int = (off / self.interval).floor();
            let next_off = (int + 1.0) * self.interval;
            let next_sample = ((next_off - off) / header.tsamp).ceil().max(1.0) as usize;
            let t1 = (t + next_sample).min(nsamp);
            let bad: Vec<bool> = if int < 0.0 || int as usize >= self.nint {
                vec![false; nchan]
            } else {
                let i = int as usize;
                map.iter().map(|ms| ms.iter().any(|&m| self.get(i, m))).collect()
            };
            if bad.iter().any(|b| *b) {
                out.push(MaskSegment { t0: t, t1, bad });
            }
            t = t1;
        }
        Ok(out)
    }

    pub fn channel_fraction(&self) -> Vec<f64> {
        (0..self.nchan()).map(|c| (0..self.nint).filter(|&i| self.get(i, c)).count() as f64 / self.nint as f64).collect()
    }
}
