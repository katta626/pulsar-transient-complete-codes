use fb_defs::Epoch;
use fb_timing::*;

fn gauss_profile(n: usize, centre: f64, width: f64) -> Vec<f64> {
    (0..n)
        .map(|i| {
            let mut d = i as f64 / n as f64 - centre;
            d -= d.round();
            (-0.5 * (d / width).powi(2)).exp()
        })
        .collect()
}

#[test]
fn fftfit_recovers_shift() {
    let t = gauss_profile(256, 0.3, 0.02);
    for shift in [0.0, 0.1234, 0.5, 0.876] {
        let mut seed = 12345u64;
        let p: Vec<f64> = gauss_profile(256, 0.3 + shift, 0.02)
            .iter()
            .map(|v| {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                3.0 * v + 0.05 * (((seed >> 33) as f64 / (1u64 << 31) as f64) - 0.5)
            })
            .collect();
        let s = fftfit(&p, &t).unwrap();
        let mut d = s.phase - (1.0 - shift).rem_euclid(1.0);
        d -= d.round();
        let alt = { let mut e = s.phase - shift; e -= e.round(); e };
        assert!(d.abs() < 1e-3 || alt.abs() < 1e-3, "shift {shift} got {}", s.phase);
        assert!(s.phase_err > 0.0 && s.phase_err < 1e-2, "err {}", s.phase_err);
        assert!((s.scale - 3.0).abs() < 0.1, "scale {}", s.scale);
    }
}

#[test]
fn polyco_round_trip() {
    let dir = std::env::temp_dir().join(format!("fbpoly_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join("polyco.dat");
    std::fs::write(
        &p,
        "1937+21    5-JAN-24  120000.00   60314.50000000000            71.0200 -0.000 -6.000\n   1234.567890  641.928000000000    7   60   3 1400.000\n 1.00000000000000000D-01 2.00000000000000000D-03-1.00000000000000000D-06\n",
    )
    .unwrap();
    let set = read_polycos(&p).unwrap();
    assert_eq!(set.len(), 1);
    let pc = &set[0];
    assert!((pc.f0 - 641.928).abs() < 1e-9);
    let t = Epoch::new(60314, 0.5 * 86400.0 + 600.0);
    let dt = 10.0;
    let expect = 1234.56789 + dt * 60.0 * 641.928 + 0.1 + 2e-3 * dt - 1e-6 * dt * dt;
    assert!((pc.phase(&t) - expect).abs() < 1e-6, "{} vs {}", pc.phase(&t), expect);
    assert!(pc.covers(&t));
    let f = pc.frequency(&t);
    assert!((f - (641.928 + (2e-3 - 2e-6 * dt) / 60.0)).abs() < 1e-12);
}

#[test]
fn tim_writer() {
    let dir = std::env::temp_dir().join(format!("fbtim_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let toas = vec![Toa { name: "obs1".into(), freq_mhz: 1400.0, epoch: Epoch::new(60000, 43200.123456789), err_us: 1.5, site: "gmrt".into(), flags: vec![] }];
    write_tim(&dir.join("x.tim"), &toas).unwrap();
    let s = std::fs::read_to_string(dir.join("x.tim")).unwrap();
    assert!(s.starts_with("FORMAT 1\nobs1 1400.000000 60000.50000142889802"), "{s}");
}
