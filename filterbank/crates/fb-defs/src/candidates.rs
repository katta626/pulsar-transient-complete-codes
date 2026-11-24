use serde::{Deserialize, Serialize};

use crate::time::Epoch;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SinglePulseCandidate {
    pub obs_id: String,
    pub dm: f64,
    pub epoch: Epoch,
    pub f_ref: f64,
    pub sample: u64,
    pub width_samples: u32,
    pub width_s: f64,
    pub snr: f64,
    pub dm_index: u32,
    pub members: u32,
    pub beam: u32,
    pub label: Option<String>,
    pub probability: Option<f64>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PeriodicCandidate {
    pub obs_id: String,
    pub f_hz: f64,
    pub fdot: f64,
    pub fddot: f64,
    pub epoch: Epoch,
    pub dm: f64,
    pub snr: f64,
    pub sigma: f64,
    pub numharm: u32,
    pub method: String,
    pub label: Option<String>,
    pub probability: Option<f64>,
}

impl PeriodicCandidate {
    pub fn period(&self) -> f64 {
        1.0 / self.f_hz
    }

    pub fn pdot(&self) -> f64 {
        -self.fdot / (self.f_hz * self.f_hz)
    }

    pub fn at_epoch(&self, t: &Epoch) -> PeriodicCandidate {
        let dt = t.seconds_since(&self.epoch);
        let mut c = self.clone();
        c.f_hz = self.f_hz + self.fdot * dt + 0.5 * self.fddot * dt * dt;
        c.fdot = self.fdot + self.fddot * dt;
        c.epoch = *t;
        c
    }
}
