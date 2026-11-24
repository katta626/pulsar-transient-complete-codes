use fb_defs::K_DM_PULSAR;
use fb_pipeline::synth::{header, Injection, SynthReader};
use fb_pipeline::{run_periodic, PeriodicConfig};

#[test]
fn finds_and_folds_dispersed_pulsar() {
    let tsamp = 1e-3;
    let n = 200_000u64;
    let h = header(64, n, tsamp, 1500.0, -4.0);
    let period = 0.333;
    let dm = 40.0;
    let injections: Vec<Injection> = (0..((n as f64 * tsamp / period) as usize))
        .map(|k| Injection { dm, t0: 1.0 + k as f64 * period, width: 6e-3, amplitude: 0.35 })
        .filter(|i| i.t0 + K_DM_PULSAR * dm * (1248f64.powi(-2) - 1500f64.powi(-2)) + 0.05 < n as f64 * tsamp)
        .collect();
    let mut r = SynthReader::new(h, &injections, &[5], 11);
    let cfg = PeriodicConfig { dm_lo: 0.0, dm_hi: 80.0, dm_step: Some(4.0), fold_top: 5, ..Default::default() };
    let res = run_periodic(&mut r, &cfg).unwrap();
    let top = &res.folded[0];
    assert!(fb_fft::is_harmonic(top.candidate.f_hz, 1.0 / period, 2e-3, 4), "top f {}", top.candidate.f_hz);
    assert!((top.candidate.dm - dm).abs() <= 12.0, "dm {}", top.candidate.dm);
    assert!(top.reduced_chi2 > 5.0, "chi2 {}", top.reduced_chi2);
    assert!(res.mask.channel_fraction()[5] > 0.5);
}
