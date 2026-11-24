use fb_defs::*;
use fb_formats::presto::{read_dat, write_dat, PrestoMask};

fn header() -> ObsHeader {
    ObsHeader {
        source_name: "TEST".into(),
        telescope: "GBT".into(),
        backend: "GUPPI".into(),
        ra_deg: Some(10.0),
        dec_deg: Some(-0.3),
        az_deg: None,
        za_deg: None,
        start: Epoch::new(59999, 12345.678901234),
        tsamp: 64e-6,
        nsamples: 10000,
        freqs: FreqAxis::uniform(1500.0, -1.0, 32),
        npol: 1,
        pol_type: PolType::Intensity,
        nbits: 8,
        scaling: Scaling::default(),
        barycentric: false,
        beam: 0,
        nbeams: 1,
        order: DataOrder::TimePolFreq,
        provenance: vec!["test".into()],
    }
}

#[test]
fn presto_mask_round_trip() {
    let h = header();
    let mut m = RfiMask::empty(&h, 0.064, "test");
    m.zap_channel(3);
    m.zap_interval(2);
    m.set(5, 10, true);
    for c in 0..32 {
        m.set(7, c, true);
    }
    let pm = PrestoMask::from_mask(&m, h.tsamp).unwrap();
    let dir = std::env::temp_dir().join(format!("fbfmt_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join("x_rfifind.mask");
    pm.write(&p).unwrap();
    let back = PrestoMask::read(&p).unwrap().to_mask_for(&h);
    let mut b = Block::new(std::sync::Arc::new(h.clone()), 0, h.nsamples as usize, vec![1.0; 32 * h.nsamples as usize]).unwrap();
    let mut b2 = b.clone();
    b.apply_mask(&m, None).unwrap();
    b2.apply_mask(&back, None).unwrap();
    assert_eq!(b.data, b2.data);
}

#[test]
fn dat_inf_round_trip() {
    let h = header();
    let s = Series { header: h.clone(), dm: 56.7, f_ref: 1500.0, data: (0..1000).map(|i| i as f32 * 0.5).collect() };
    let dir = std::env::temp_dir().join(format!("fbdat_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    write_dat(&dir.join("obs_DM56.70"), &s).unwrap();
    let back = read_dat(&dir.join("obs_DM56.70.dat")).unwrap();
    assert_eq!(back.data, s.data);
    assert!((back.dm - 56.7).abs() < 1e-9);
    assert!(back.header.start.seconds_since(&h.start).abs() < 1e-6);
    assert!((back.header.dec_deg.unwrap() + 0.3).abs() < 1e-6);
    assert!((back.header.tsamp - 64e-6).abs() < 1e-18);
}
