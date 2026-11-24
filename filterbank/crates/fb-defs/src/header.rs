use serde::{Deserialize, Serialize};

use crate::error::{FbError, Result};
use crate::time::Epoch;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FreqAxis {
    pub mhz: Vec<f64>,
    pub chan_bw: f64,
}

impl FreqAxis {
    pub fn uniform(fch1: f64, foff: f64, nchan: usize) -> Self {
        Self { mhz: (0..nchan).map(|i| fch1 + foff * i as f64).collect(), chan_bw: foff }
    }

    pub fn from_centres(mhz: Vec<f64>, chan_bw: Option<f64>) -> Result<Self> {
        if mhz.is_empty() {
            return Err(FbError::Invalid("frequency axis has no channels".into()));
        }
        let bw = match chan_bw {
            Some(b) if b != 0.0 => b,
            _ if mhz.len() > 1 => median_step(&mhz),
            _ => return Err(FbError::Invalid("single-channel axis needs a channel width".into())),
        };
        Ok(Self { mhz, chan_bw: bw })
    }

    pub fn nchan(&self) -> usize {
        self.mhz.len()
    }

    pub fn fch1(&self) -> f64 {
        self.mhz[0]
    }

    pub fn foff(&self) -> f64 {
        self.chan_bw
    }

    pub fn is_uniform(&self, rel_tol: f64) -> bool {
        let tol = self.chan_bw.abs() * rel_tol;
        self.mhz.iter().enumerate().all(|(i, f)| (f - (self.mhz[0] + self.chan_bw * i as f64)).abs() <= tol)
    }

    pub fn is_descending(&self) -> bool {
        self.chan_bw < 0.0
    }

    pub fn lo(&self) -> f64 {
        self.mhz.iter().cloned().fold(f64::INFINITY, f64::min) - self.chan_bw.abs() / 2.0
    }

    pub fn hi(&self) -> f64 {
        self.mhz.iter().cloned().fold(f64::NEG_INFINITY, f64::max) + self.chan_bw.abs() / 2.0
    }

    pub fn centre(&self) -> f64 {
        0.5 * (self.lo() + self.hi())
    }

    pub fn bandwidth(&self) -> f64 {
        self.hi() - self.lo()
    }

    pub fn edges(&self, chan: usize) -> (f64, f64) {
        let h = self.chan_bw.abs() / 2.0;
        (self.mhz[chan] - h, self.mhz[chan] + h)
    }

    pub fn reversed(&self) -> Self {
        let mut mhz = self.mhz.clone();
        mhz.reverse();
        Self { mhz, chan_bw: -self.chan_bw }
    }

    pub fn same_as(&self, other: &FreqAxis, rel_tol: f64) -> bool {
        let tol = self.chan_bw.abs() * rel_tol;
        self.nchan() == other.nchan()
            && self.mhz.iter().zip(&other.mhz).all(|(a, b)| (a - b).abs() <= tol)
    }
}

fn median_step(mhz: &[f64]) -> f64 {
    let mut d: Vec<f64> = mhz.windows(2).map(|w| w[1] - w[0]).collect();
    d.sort_by(|a, b| a.partial_cmp(b).unwrap());
    d[d.len() / 2]
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PolType {
    Intensity,
    AABB,
    AABBCRCI,
    IQUV,
    Unknown,
}

impl PolType {
    pub fn parse(s: &str, npol: usize) -> Self {
        match s.trim().to_ascii_uppercase().as_str() {
            "AA+BB" | "INTEN" | "I" => PolType::Intensity,
            "AABB" => PolType::AABB,
            "AABBCRCI" => PolType::AABBCRCI,
            "IQUV" => PolType::IQUV,
            _ if npol == 1 => PolType::Intensity,
            _ => PolType::Unknown,
        }
    }

    pub fn label(&self) -> &'static str {
        match self {
            PolType::Intensity => "AA+BB",
            PolType::AABB => "AABB",
            PolType::AABBCRCI => "AABBCRCI",
            PolType::IQUV => "IQUV",
            PolType::Unknown => "UNKNOWN",
        }
    }

    pub fn total_intensity_terms(&self, npol: usize) -> Vec<usize> {
        match self {
            PolType::AABB | PolType::AABBCRCI if npol >= 2 => vec![0, 1],
            _ => vec![0],
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DataOrder {
    TimePolFreq,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Scaling {
    pub zero_off: f64,
    pub signed: bool,
    pub applied: bool,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ObsHeader {
    pub source_name: String,
    pub telescope: String,
    pub backend: String,
    pub ra_deg: Option<f64>,
    pub dec_deg: Option<f64>,
    pub az_deg: Option<f64>,
    pub za_deg: Option<f64>,
    pub start: Epoch,
    pub tsamp: f64,
    pub nsamples: u64,
    pub freqs: FreqAxis,
    pub npol: usize,
    pub pol_type: PolType,
    pub nbits: u32,
    pub scaling: Scaling,
    pub barycentric: bool,
    pub beam: u32,
    pub nbeams: u32,
    pub order: DataOrder,
    pub provenance: Vec<String>,
}

impl ObsHeader {
    pub fn nchan(&self) -> usize {
        self.freqs.nchan()
    }

    pub fn duration(&self) -> f64 {
        self.nsamples as f64 * self.tsamp
    }

    pub fn epoch_of_sample(&self, sample: u64) -> Epoch {
        self.start.add_seconds(sample as f64 * self.tsamp)
    }

    pub fn sample_of_epoch(&self, t: &Epoch) -> f64 {
        t.seconds_since(&self.start) / self.tsamp
    }

    pub fn observation_id(&self) -> String {
        format!(
            "{}|{}|{}|{:.3}|beam{}",
            self.telescope.trim(),
            self.source_name.trim(),
            self.start.imjd,
            self.start.sec,
            self.beam
        )
    }

    pub fn with_step(&self, step: impl Into<String>) -> Self {
        let mut h = self.clone();
        h.provenance.push(step.into());
        h
    }

    pub fn check_same_observation(&self, other: &ObsHeader) -> Result<()> {
        if self.observation_id() != other.observation_id() {
            return Err(FbError::Incompatible(format!(
                "different observations: '{}' vs '{}'",
                self.observation_id(),
                other.observation_id()
            )));
        }
        Ok(())
    }

    pub fn check_same_freqs(&self, other: &FreqAxis) -> Result<()> {
        if !self.freqs.same_as(other, 1e-3) {
            return Err(FbError::Incompatible(format!(
                "frequency axes differ: {} ch {:.4}..{:.4} MHz vs {} ch {:.4}..{:.4} MHz",
                self.nchan(),
                self.freqs.lo(),
                self.freqs.hi(),
                other.nchan(),
                other.lo(),
                other.hi()
            )));
        }
        Ok(())
    }

    pub fn ra_hms(&self) -> Option<(i32, i32, f64)> {
        self.ra_deg.map(|d| {
            let h = d.rem_euclid(360.0) / 15.0;
            let hh = h.floor();
            let m = (h - hh) * 60.0;
            let mm = m.floor();
            (hh as i32, mm as i32, (m - mm) * 60.0)
        })
    }

    pub fn dec_dms(&self) -> Option<(bool, i32, i32, f64)> {
        self.dec_deg.map(|d| {
            let neg = d < 0.0;
            let a = d.abs();
            let dd = a.floor();
            let m = (a - dd) * 60.0;
            let mm = m.floor();
            (neg, dd as i32, mm as i32, (m - mm) * 60.0)
        })
    }
}

pub fn parse_sexagesimal(text: &str) -> Option<(bool, f64)> {
    let t = text.trim();
    if t.is_empty() {
        return None;
    }
    let neg = t.starts_with('-');
    let body = t.trim_start_matches(['+', '-']);
    let parts: Vec<f64> = body.split(':').map(|p| p.trim().parse::<f64>()).collect::<std::result::Result<_, _>>().ok()?;
    let mut v = 0.0;
    for (i, p) in parts.iter().enumerate() {
        v += p / 60f64.powi(i as i32);
    }
    Some((neg, v))
}

pub fn ra_from_hms(text: &str) -> Option<f64> {
    parse_sexagesimal(text).map(|(_, h)| h * 15.0)
}

pub fn dec_from_dms(text: &str) -> Option<f64> {
    parse_sexagesimal(text).map(|(neg, d)| if neg { -d } else { d })
}
