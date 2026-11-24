use std::sync::Arc;

use crate::error::{FbError, Result};
use crate::header::ObsHeader;
use crate::mask::RfiMask;
use crate::time::Epoch;

#[derive(Debug, Clone)]
pub struct Block {
    pub header: Arc<ObsHeader>,
    pub start_sample: u64,
    pub nsamp: usize,
    pub data: Vec<f32>,
}

impl Block {
    pub fn new(header: Arc<ObsHeader>, start_sample: u64, nsamp: usize, data: Vec<f32>) -> Result<Self> {
        let expect = nsamp * header.npol * header.nchan();
        if data.len() != expect {
            return Err(FbError::Invalid(format!("block has {} values, expected {}", data.len(), expect)));
        }
        Ok(Self { header, start_sample, nsamp, data })
    }

    pub fn nchan(&self) -> usize {
        self.header.nchan()
    }

    pub fn npol(&self) -> usize {
        self.header.npol
    }

    pub fn start_epoch(&self) -> Epoch {
        self.header.epoch_of_sample(self.start_sample)
    }

    pub fn index(&self, t: usize, p: usize, c: usize) -> usize {
        (t * self.npol() + p) * self.nchan() + c
    }

    pub fn spectrum(&self, t: usize, p: usize) -> &[f32] {
        let i = self.index(t, p, 0);
        &self.data[i..i + self.nchan()]
    }

    pub fn total_intensity(&self) -> Block {
        let nchan = self.nchan();
        let terms = self.header.pol_type.total_intensity_terms(self.npol());
        if self.npol() == 1 {
            return self.clone();
        }
        let mut out = vec![0f32; self.nsamp * nchan];
        for t in 0..self.nsamp {
            for &p in &terms {
                let s = self.spectrum(t, p);
                for (o, v) in out[t * nchan..(t + 1) * nchan].iter_mut().zip(s) {
                    *o += *v;
                }
            }
        }
        let mut h = (*self.header).clone();
        h.npol = 1;
        h.pol_type = crate::header::PolType::Intensity;
        h.provenance.push("total_intensity".into());
        Block { header: Arc::new(h), start_sample: self.start_sample, nsamp: self.nsamp, data: out }
    }

    pub fn apply_mask(&mut self, mask: &RfiMask, fill: Option<&[f32]>) -> Result<usize> {
        mask.check_applies_to(&self.header)?;
        let nchan = self.nchan();
        let npol = self.npol();
        let segments = mask.to_indices(&self.header, self.start_sample, self.nsamp)?;
        let mut zapped = 0usize;
        for seg in &segments {
            for t in seg.t0..seg.t1 {
                for p in 0..npol {
                    let base = (t * npol + p) * nchan;
                    for c in 0..nchan {
                        if seg.bad[c] {
                            self.data[base + c] = fill.map(|f| f[c]).unwrap_or(0.0);
                            zapped += 1;
                        }
                    }
                }
            }
        }
        Ok(zapped)
    }
}
