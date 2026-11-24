use serde::{Deserialize, Serialize};

use crate::consts::SECONDS_PER_DAY;

#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
pub struct Epoch {
    pub imjd: i64,
    pub sec: f64,
}

impl Epoch {
    pub fn new(imjd: i64, sec: f64) -> Self {
        Self { imjd, sec }.normalized()
    }

    pub fn from_mjd(mjd: f64) -> Self {
        let imjd = mjd.floor() as i64;
        Self::new(imjd, (mjd - imjd as f64) * SECONDS_PER_DAY)
    }

    pub fn from_parts(imjd: i64, smjd: f64, offs: f64) -> Self {
        Self::new(imjd, smjd + offs)
    }

    pub fn normalized(self) -> Self {
        let days = (self.sec / SECONDS_PER_DAY).floor();
        Self { imjd: self.imjd + days as i64, sec: self.sec - days * SECONDS_PER_DAY }
    }

    pub fn mjd(&self) -> f64 {
        self.imjd as f64 + self.sec / SECONDS_PER_DAY
    }

    pub fn add_seconds(&self, dt: f64) -> Self {
        Self::new(self.imjd, self.sec + dt)
    }

    pub fn seconds_since(&self, other: &Epoch) -> f64 {
        (self.imjd - other.imjd) as f64 * SECONDS_PER_DAY + (self.sec - other.sec)
    }
}
