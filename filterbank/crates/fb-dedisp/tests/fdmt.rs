use std::sync::Arc;
use std::time::Instant;

use fb_defs::*;

fn synth(nchan: usize, nt: usize, dm: f64, t0: f64, tsamp: f64) -> Block {
    let h = ObsHeader {
        source_name: "S".into(),
        telescope: "Fake".into(),
        backend: "FAKE".into(),
        ra_deg: None,
        dec_deg: None,
        az_deg: None,
        za_deg: None,
        start: Epoch::new(60000, 0.0),
        tsamp,
        nsamples: nt as u64,
        freqs: FreqAxis::uniform(1500.0, -300.0 / nchan as f64, nchan),
        npol: 1,
        pol_type: PolType::Intensity,
        nbits: 32,
        scaling: Scaling::default(),
        barycentric: false,
        beam: 0,
        nbeams: 1,
        order: DataOrder::TimePolFreq,
        provenance: vec![],
    };
    let mut seed = 7u64;
    let mut data = vec![0f32; nt * nchan];
    for v in data.iter_mut() {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *v = ((seed >> 40) as f32 / (1u64 << 24) as f32 - 0.5) * 3.46;
    }
    let f_ref = h.freqs.fch1();
    for c in 0..nchan {
        let s = ((t0 + K_DM_PULSAR * dm * (h.freqs.mhz[c].powi(-2) - f_ref.powi(-2))) / tsamp).round() as usize;
        for k in 0..4 {
            if s + k < nt {
                data[(s + k) * nchan + c] += 2.0;
            }
        }
    }
    Block::new(Arc::new(h), 0, nt, data).unwrap()
}

#[test]
fn fdmt_matches_brute_force() {
    let tsamp = 256e-6;
    let b = synth(256, 1 << 15, 300.0, 1.0, tsamp);
    let plan = DmPlan::linear(b.header.freqs.clone(), tsamp, 0.0, 500.0, 0.5).unwrap();
    let t = Instant::now();
    let bf = fb_dedisp::brute_force(&b, &plan, None).unwrap();
    let t_bf = t.elapsed();
    let t = Instant::now();
    let fd = fb_dedisp::fdmt(&b, &plan, None).unwrap();
    let t_fd = t.elapsed();
    assert_eq!(bf.nout, fd.nout);
    let peak = |o: &fb_dedisp::DedispOutput| {
        let (i, v) = o.data.iter().enumerate().fold((0, f32::MIN), |a, (i, v)| if *v > a.1 { (i, *v) } else { a });
        (i / o.nout, i % o.nout, v)
    };
    let (d1, t1, v1) = peak(&bf);
    let (d2, t2, v2) = peak(&fd);
    println!("brute {:?} dm {} t {} v {v1};  fdmt {:?} dm {} t {} v {v2}", t_bf, plan.dms[d1], t1, t_fd, plan.dms[d2], t2);
    assert!((plan.dms[d1] - 300.0).abs() < 3.0 && (plan.dms[d2] - 300.0).abs() < 3.0);
    assert!((t1 as i64 - t2 as i64).abs() <= 3, "time {t1} vs {t2}");
    assert!(v2 > 0.85 * v1, "fdmt peak {v2} vs brute {v1}");
}

#[test]
fn chunked_fdmt_is_identical() {
    let tsamp = 256e-6;
    let b = synth(64, 1 << 14, 200.0, 1.0, tsamp);
    let plan = DmPlan::linear(b.header.freqs.clone(), tsamp, 0.0, 300.0, 1.0).unwrap();
    let full = fb_dedisp::fdmt::fdmt_chunked(&b, &plan, None, usize::MAX / 4).unwrap();
    let small = fb_dedisp::fdmt::fdmt_chunked(&b, &plan, None, 1 << 18).unwrap();
    assert_eq!(full.nout, small.nout);
    let maxdiff = full.data.iter().zip(&small.data).map(|(a, b)| (a - b).abs()).fold(0f32, f32::max);
    assert!(maxdiff < 1e-3, "max diff {maxdiff}");
}

#[test]
#[ignore]
fn bench_4096_channels() {
    let tsamp = 64e-6;
    let b = synth(4096, 1 << 15, 200.0, 0.5, tsamp);
    let plan = DmPlan::linear(b.header.freqs.clone(), tsamp, 0.0, 250.0, 0.25).unwrap();
    let t = Instant::now();
    let bf = fb_dedisp::brute_force(&b, &plan, None).unwrap();
    let t_bf = t.elapsed();
    let t = Instant::now();
    let sb = fb_dedisp::subband(&b, &plan, 64, None).unwrap();
    let t_sb = t.elapsed();
    let t = Instant::now();
    let fd = fb_dedisp::fdmt(&b, &plan, None).unwrap();
    let t_fd = t.elapsed();
    println!("4096 ch, {} DMs, {} samples ({:.2} s of data): brute {:?}, subband {:?}, fdmt {:?}", plan.ndm(), 1 << 15, (1 << 15) as f64 * tsamp, t_bf, t_sb, t_fd);
    let _ = (bf.nout, sb.nout, fd.nout);
}
