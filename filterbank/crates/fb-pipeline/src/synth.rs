use std::sync::Arc;

use fb_defs::{Block, DataOrder, Epoch, FreqAxis, ObsHeader, PolType, Result, Scaling, K_DM_PULSAR};
use fb_io::Gap;
use fb_io::Reader;

pub struct SynthReader {
    header: Arc<ObsHeader>,
    data: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct Injection {
    pub dm: f64,
    pub t0: f64,
    pub width: f64,
    pub amplitude: f32,
}

pub struct Lcg(pub u64);

impl Lcg {
    pub fn next_f32(&mut self) -> f32 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 40) as f32) / (1u64 << 24) as f32
    }

    pub fn gauss(&mut self) -> f32 {
        let u1 = self.next_f32().max(1e-7);
        let u2 = self.next_f32();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f32::consts::PI * u2).cos()
    }
}

pub fn header(nchan: usize, nsamp: u64, tsamp: f64, fch1: f64, foff: f64) -> ObsHeader {
    ObsHeader {
        source_name: "SYNTH".into(),
        telescope: "Fake".into(),
        backend: "FAKE".into(),
        ra_deg: Some(0.0),
        dec_deg: Some(0.0),
        az_deg: None,
        za_deg: None,
        start: Epoch::new(60000, 0.0),
        tsamp,
        nsamples: nsamp,
        freqs: FreqAxis::uniform(fch1, foff, nchan),
        npol: 1,
        pol_type: PolType::Intensity,
        nbits: 32,
        scaling: Scaling::default(),
        barycentric: false,
        beam: 0,
        nbeams: 1,
        order: DataOrder::TimePolFreq,
        provenance: vec!["synthetic".into()],
    }
}

impl SynthReader {
    pub fn new(h: ObsHeader, injections: &[Injection], rfi_channels: &[usize], seed: u64) -> Self {
        let nchan = h.nchan();
        let n = h.nsamples as usize;
        let mut rng = Lcg(seed);
        let mut data: Vec<f32> = (0..n * nchan).map(|_| 10.0 + rng.gauss()).collect();
        let f_ref = h.freqs.mhz.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        for inj in injections {
            for c in 0..nchan {
                let delay = K_DM_PULSAR * inj.dm * (h.freqs.mhz[c].powi(-2) - f_ref.powi(-2));
                let tc = inj.t0 + delay;
                let s0 = ((tc - 3.0 * inj.width) / h.tsamp).floor().max(0.0) as usize;
                let s1 = (((tc + 3.0 * inj.width) / h.tsamp).ceil() as usize).min(n);
                for s in s0..s1 {
                    let t = s as f64 * h.tsamp;
                    let g = (-0.5 * ((t - tc) / inj.width).powi(2)).exp() as f32;
                    data[s * nchan + c] += inj.amplitude * g;
                }
            }
        }
        for &c in rfi_channels {
            for s in 0..n {
                data[s * nchan + c] += 30.0 * rng.gauss().abs() + 20.0 * ((s as f32) * 0.01).sin();
            }
        }
        Self { header: Arc::new(h), data }
    }

    pub fn data(&self) -> &[f32] {
        &self.data
    }
}

impl Reader for SynthReader {
    fn header(&self) -> Arc<ObsHeader> {
        self.header.clone()
    }

    fn gaps(&self) -> &[Gap] {
        &[]
    }

    fn read(&mut self, start: u64, nsamp: usize) -> Result<Block> {
        let nchan = self.header.nchan();
        let end = (start + nsamp as u64).min(self.header.nsamples);
        let n = end.saturating_sub(start) as usize;
        let d = self.data[start as usize * nchan..end as usize * nchan].to_vec();
        Block::new(self.header.clone(), start, n, d)
    }
}
