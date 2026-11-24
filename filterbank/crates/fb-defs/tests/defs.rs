use std::sync::Arc;

use fb_defs::*;

fn header(nchan: usize, fch1: f64, foff: f64) -> ObsHeader {
    ObsHeader {
        source_name: "TEST".into(),
        telescope: "GMRT".into(),
        backend: "FAKE".into(),
        ra_deg: Some(83.63),
        dec_deg: Some(-0.25),
        az_deg: None,
        za_deg: None,
        start: Epoch::new(60000, 3600.25),
        tsamp: 1e-3,
        nsamples: 1000,
        freqs: FreqAxis::uniform(fch1, foff, nchan),
        npol: 1,
        pol_type: PolType::Intensity,
        nbits: 8,
        scaling: Scaling::default(),
        barycentric: false,
        beam: 1,
        nbeams: 1,
        order: DataOrder::TimePolFreq,
        provenance: vec![],
    }
}

#[test]
fn epoch_keeps_precision() {
    let e = Epoch::new(60000, 86399.999_999_9);
    let f = e.add_seconds(1e-7);
    assert_eq!(f.imjd, 60001);
    assert!((f.seconds_since(&e) - 1e-7).abs() < 1e-12);
}

#[test]
fn mask_follows_band_flip() {
    let h = header(8, 1500.0, -25.0);
    let mut m = RfiMask::empty(&h, 0.1, "test");
    m.zap_freq_range(1440.0, 1460.0);
    let flipped = ObsHeader { freqs: h.freqs.reversed(), ..h.clone() };
    let mut b = Block::new(Arc::new(flipped.clone()), 0, 10, vec![1.0; 80]).unwrap();
    b.apply_mask(&m, None).unwrap();
    let bad: Vec<usize> = (0..8).filter(|&c| b.data[c] == 0.0).collect();
    let freqs: Vec<f64> = bad.iter().map(|&c| flipped.freqs.mhz[c]).collect();
    assert_eq!(freqs, vec![1450.0]);
}

#[test]
fn mask_rejects_other_observation() {
    let h = header(8, 1500.0, -25.0);
    let m = RfiMask::empty(&h, 0.1, "test");
    let other = ObsHeader { source_name: "OTHER".into(), ..h };
    let mut b = Block::new(Arc::new(other), 0, 10, vec![1.0; 80]).unwrap();
    assert!(b.apply_mask(&m, None).is_err());
}

#[test]
fn dm_delays_shared_table() {
    let h = header(4, 1500.0, -100.0);
    let p = DmPlan::linear(h.freqs.clone(), 1e-3, 0.0, 100.0, 50.0).unwrap();
    let t = p.delay_table();
    assert_eq!(t.len(), 12);
    assert_eq!(t[0], 0);
    let expect = (K_DM_PULSAR * 100.0 * (1200f64.powi(-2) - 1500f64.powi(-2)) / 1e-3).round() as i32;
    assert_eq!(t[11], expect);
}

#[test]
fn dm_tolerance_plan_increases() {
    let h = header(256, 1500.0, -1.0);
    let p = DmPlan::tolerance(h.freqs, 64e-6, 0.0, 500.0, 1e-3, 1.25).unwrap();
    assert!(p.ndm() > 10);
    assert!(p.dms.windows(2).all(|w| w[1] > w[0]));
}
