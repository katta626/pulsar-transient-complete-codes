use fb_defs::{Epoch, Series};
use fb_fft::fdas::{self, FdasParams};
use fb_fft::FftSearchParams;
use fb_pipeline::synth::{header, Lcg};

fn pulsar(n: usize, tsamp: f64, f0: f64, fdot: f64, fddot: f64, duty: f64, amp: f32, seed: u64) -> Series {
    let mut h = header(1, n as u64, tsamp, 1400.0, -1.0);
    h.start = Epoch::new(60000, 0.0);
    let mut rng = Lcg(seed);
    let t_mid = n as f64 * tsamp / 2.0;
    let data: Vec<f32> = (0..n)
        .map(|i| {
            let tt = i as f64 * tsamp - t_mid;
            let ph = (f0 * tt + fdot / 2.0 * tt * tt + fddot / 6.0 * tt * tt * tt).rem_euclid(1.0);
            rng.gauss() + if (ph - 0.5).abs() < duty / 2.0 { amp } else { 0.0 }
        })
        .collect();
    Series { header: h, dm: 30.0, f_ref: 1400.0, data }
}

#[test]
fn response_matches_sinc() {
    let k = fdas::response(0.0, 0.0, 0.0, 16);
    assert!((k[16].norm() - 1.0).abs() < 1e-4);
    assert!(k.iter().enumerate().all(|(i, v)| i == 16 || v.norm() < 1e-4));
    let h = fdas::response(0.0, 0.0, 0.5, 16);
    assert!((h[16].norm() - 2.0 / std::f32::consts::PI).abs() < 1e-3);
    let z = fdas::response(40.0, 0.0, 0.0, 40);
    let p: f32 = z.iter().map(|v| v.norm_sqr()).sum();
    assert!((p - 1.0).abs() < 0.02, "Parseval {p}");
}

#[test]
fn fdas_recovers_accelerated_pulsar() {
    let n = 1 << 18;
    let ts = 1e-3;
    let t = n as f64 * ts;
    let (f0, z0) = (7.3, 24.0);
    let s = pulsar(n, ts, f0, z0 / (t * t), 0.0, 0.08, 0.1, 3);
    let params = FftSearchParams { f_min: 5.0, f_max: 10.0, max_harm: 8, ..Default::default() };
    let flat = fdas::search_series(&s, &[], &FdasParams { zmax: 0.0, refine: 0, ..Default::default() }, &params);
    let fp = FdasParams { zmax: 120.0, refine: 5, ..Default::default() };
    let c = fdas::search_series(&s, &[], &fp, &params);
    let top = &c[0];
    assert!((top.f_hz - f0).abs() < 1.0 / t, "f {} vs {}", top.f_hz, f0);
    assert!((top.fdot * t * t - z0).abs() < 1.0, "z {}", top.fdot * t * t);
    let flat_best = flat.first().map(|x| x.sigma).unwrap_or(0.0);
    assert!(top.sigma > flat_best + 3.0, "accel {} vs flat {}", top.sigma, flat_best);
}

#[test]
fn strong_drift_is_not_whitened_away() {
    let n = 1 << 19;
    let ts = 64e-6;
    let t = n as f64 * ts;
    let (f0, z0) = (23.4567, 180.0);
    let s = pulsar(n, ts, f0, z0 / (t * t), 0.0, 0.1, 0.12, 11);
    let (mut re, mut im) = (0f64, 0f64);
    for (i, v) in s.data.iter().enumerate() {
        let tt = i as f64 * ts - t / 2.0;
        let a = std::f64::consts::TAU * (f0 * tt + 0.5 * z0 / (t * t) * tt * tt);
        re += *v as f64 * a.cos();
        im -= *v as f64 * a.sin();
    }
    let ideal = (re * re + im * im) / n as f64;
    let spec = fb_fft::prepare_width(&s, &[], 8 * 200);
    let got = fdas::power_at(&spec, f0 * t, z0, 0.0) as f64;
    assert!(got > 0.8 * ideal, "power {got} of ideal {ideal}");
}

#[test]
fn jerk_search_recovers_w() {
    let n = 1 << 18;
    let ts = 1e-3;
    let t = n as f64 * ts;
    let (f0, z0, w0) = (7.3, 10.0, 60.0);
    let s = pulsar(n, ts, f0, z0 / (t * t), w0 / (t * t * t), 0.1, 0.1, 5);
    let params = FftSearchParams { f_min: 6.0, f_max: 9.0, max_harm: 4, ..Default::default() };
    let accel_only = fdas::search_series(&s, &[], &FdasParams { zmax: 100.0, refine: 0, ..Default::default() }, &params);
    let fp = FdasParams { zmax: 100.0, wmax: 300.0, dw: 20.0, refine: 5, ..Default::default() };
    let c = fdas::search_series(&s, &[], &fp, &params);
    let top = &c[0];
    assert!((top.f_hz - f0).abs() < 1.0 / t, "f {}", top.f_hz);
    assert!((top.fddot * t * t * t - w0).abs() < 15.0, "w {}", top.fddot * t * t * t);
    let acc = accel_only.first().map(|x| x.sigma).unwrap_or(0.0);
    assert!(top.sigma > acc, "jerk {} vs accel {}", top.sigma, acc);
}

#[test]
#[ignore]
fn bench_fdas_full_range() {
    let n = 1 << 20;
    let ts = 5e-4;
    let t = n as f64 * ts;
    let s = pulsar(n, ts, 31.7, 50.0 / (t * t), 0.0, 0.05, 0.1, 9);
    let params = FftSearchParams { f_min: 0.5, f_max: 500.0, max_harm: 16, ..Default::default() };
    let t0 = std::time::Instant::now();
    let c = fdas::search_series(&s, &[], &FdasParams { zmax: 200.0, refine: 20, ..Default::default() }, &params);
    println!("fdas zmax=200 harm=16 n=2^20: {:.2} s, top f={:.5} sigma={:.1}", t0.elapsed().as_secs_f64(), c[0].f_hz, c[0].sigma);
}
