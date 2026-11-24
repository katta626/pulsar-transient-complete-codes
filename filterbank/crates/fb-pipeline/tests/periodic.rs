use fb_defs::{Epoch, Series, SPEED_OF_LIGHT};
use fb_pipeline::synth::{header, Lcg};

fn pulsar(n: usize, tsamp: f64, period: f64, duty: f64, amp: f32, accel: f64, seed: u64) -> Series {
    let mut h = header(1, n as u64, tsamp, 1400.0, -1.0);
    h.start = Epoch::new(60000, 0.0);
    let mut rng = Lcg(seed);
    let t_mid = n as f64 * tsamp / 2.0;
    let data: Vec<f32> = (0..n)
        .map(|i| {
            let t = i as f64 * tsamp;
            let tt = t - t_mid;
            let delay = accel * tt * tt / (2.0 * SPEED_OF_LIGHT);
            let ph = ((t - delay) / period).rem_euclid(1.0);
            let d = (ph - 0.5).abs();
            rng.gauss() + if d < duty / 2.0 { amp } else { 0.0 }
        })
        .collect();
    Series { header: h, dm: 10.0, f_ref: 1400.0, data }
}

#[test]
fn fft_finds_pulsar_with_harmonics() {
    let s = pulsar(1 << 20, 1e-3, 0.25, 0.04, 0.15, 0.0, 1);
    let params = fb_fft::FftSearchParams { f_min: 0.5, f_max: 100.0, ..Default::default() };
    let c = fb_fft::search_series(&s, &[0.0], &[], &params);
    assert!(!c.is_empty());
    assert!(fb_fft::is_harmonic(c[0].f_hz, 4.0, 1e-3, 8), "top {} Hz", c[0].f_hz);
    assert!(c[0].sigma > 8.0, "sigma {}", c[0].sigma);
    let sifted = fb_fft::sift(&c, s.len() as f64 * 1e-3, 1, 8);
    assert!((sifted[0].f_hz - 4.0).abs() < 0.01, "sifted top {}", sifted[0].f_hz);
    assert!(sifted.iter().filter(|x| x.sigma > 8.0).count() <= 2, "harmonics removed: {}", sifted.len());
}

#[test]
fn accel_search_recovers_binary() {
    let a = 60.0;
    let s = pulsar(1 << 19, 5e-4, 0.005, 0.1, 0.25, a, 7);
    let t_obs = s.len() as f64 * 5e-4;
    let params = fb_fft::FftSearchParams { f_min: 50.0, f_max: 400.0, max_harm: 8, ..Default::default() };
    let none = fb_fft::search_series(&s, &[0.0], &[], &params);
    let trials = fb_fft::accel_trials(t_obs, 5e-4, 100.0, 5e-4, 1.1);
    let acc = fb_fft::search_series(&s, &trials, &[], &params);
    let best = acc.iter().find(|c| (c.f_hz - 200.0).abs() < 0.1).expect("found at 200 Hz");
    let best_none = none.iter().find(|c| (c.f_hz - 200.0).abs() < 0.1).map(|c| c.sigma).unwrap_or(0.0);
    assert!(best.sigma > best_none + 3.0, "accel {} vs none {}", best.sigma, best_none);
    let est = -best.fdot * SPEED_OF_LIGHT / best.f_hz;
    assert!((est - a).abs() < 25.0, "accel estimate {est}");
}

#[test]
fn ffa_finds_long_period() {
    let s = pulsar(600_000, 1e-3, 4.2, 0.02, 0.5, 0.0, 3);
    let params = fb_ffa::FfaParams { period_min: 2.0, period_max: 8.0, bins_min: 200, bins_max: 220, snr_threshold: 8.0, duty_max: 0.1 };
    let c = fb_ffa::search(&s, &params);
    assert!(!c.is_empty());
    let p = c[0].period();
    let ok = [4.2, 8.4, 2.1].iter().any(|x| (p - x).abs() < 0.01);
    assert!(ok && c.iter().any(|x| (x.period() - 4.2).abs() < 0.005), "top period {p}");
    assert!(c[0].snr > 15.0, "snr {}", c[0].snr);
}

#[test]
fn fold_and_optimise() {
    let s = pulsar(400_000, 1e-3, 0.5, 0.05, 0.2, 0.0, 9);
    let fp = fb_fold::FoldParams { f: 2.0 + 5e-4, fdot: 0.0, fddot: 0.0, epoch: s.header.start, dm: 10.0, nbins: 64 };
    let good = fb_fold::fold_series(&s, &fb_fold::FoldParams { f: 2.0, ..fp }, 8).reduced_chi2(1.0);
    let bad = fb_fold::fold_series(&s, &fb_fold::FoldParams { f: 2.013, ..fp }, 8).reduced_chi2(1.0);
    assert!(good > 10.0 * bad.max(1.0), "chi2 {good} vs {bad}");
    let (best, chi) = fb_fold::optimise(&s, &fp, 1e-3, 81, 1, 0.0, 1.0);
    assert!((best.f - 2.0).abs() < 5e-5 && chi > 0.9 * good, "opt f {} chi {chi}", best.f);
}
