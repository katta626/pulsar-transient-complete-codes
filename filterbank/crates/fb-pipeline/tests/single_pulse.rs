use fb_pipeline::synth::{header, Injection, SynthReader};
use fb_pipeline::{run_single_pulse, Dedisperser, SinglePulseConfig};

fn run(ded: Dedisperser) -> Vec<fb_defs::SinglePulseCandidate> {
    let h = header(128, 120_000, 256e-6, 1500.0, -2.0);
    let inj = vec![
        Injection { dm: 150.0, t0: 7.0, width: 2e-3, amplitude: 1.2 },
        Injection { dm: 400.0, t0: 20.0, width: 4e-3, amplitude: 1.0 },
    ];
    let mut r = SynthReader::new(h, &inj, &[17, 90], 42);
    let mut cfg = SinglePulseConfig { dm_lo: 0.0, dm_hi: 600.0, dedisperser: ded, gulp: 40_000, ..Default::default() };
    cfg.search.threshold = 8.0;
    cfg.search.max_width = 256;
    let res = run_single_pulse(&mut r, &cfg).unwrap();
    assert!(res.mask.channel_fraction()[17] > 0.5, "RFI channel 17 masked");
    assert!(res.mask.channel_fraction()[90] > 0.5, "RFI channel 90 masked");
    res.candidates
}

#[test]
fn recovers_injected_pulses() {
    for ded in [Dedisperser::BruteForce, Dedisperser::Subband(16)] {
        let c = run(ded);
        let best = |t: f64| {
            c.iter()
                .filter(|x| (x.epoch.seconds_since(&fb_defs::Epoch::new(60000, 0.0)) - t).abs() < 0.05)
                .max_by(|a, b| a.snr.partial_cmp(&b.snr).unwrap())
                .cloned()
        };
        let a = best(7.0).expect("pulse at 7 s found");
        let b = best(20.0).expect("pulse at 20 s found");
        assert!((a.dm - 150.0).abs() < 8.0, "{ded:?} dm {}", a.dm);
        assert!((b.dm - 400.0).abs() < 15.0, "{ded:?} dm {}", b.dm);
        assert!(a.snr > 15.0 && b.snr > 15.0, "{ded:?} snr {} {}", a.snr, b.snr);
        let false_alarms = c.iter().filter(|x| x.snr > 10.0 && best(7.0).map(|a| (a.sample as i64 - x.sample as i64).abs() > 2000).unwrap_or(true) && best(20.0).map(|b| (b.sample as i64 - x.sample as i64).abs() > 2000).unwrap_or(true)).count();
        assert!(false_alarms == 0, "{ded:?} false alarms above 10 sigma: {false_alarms}");
    }
}
