use fb_defs::*;
use fb_io::bits::{pack, unpack, BitOrder};
use fb_io::{open, Reader, SigprocWriter, WriteOptions};

#[test]
fn bit_round_trip_both_orders() {
    for nbits in [1u32, 2, 4, 8] {
        let max = ((1u32 << nbits) - 1) as f32;
        let vals: Vec<f32> = (0..64).map(|i| (i as f32) % (max + 1.0)).collect();
        for order in [BitOrder::MsbFirst, BitOrder::LsbFirst] {
            let p = pack(&vals, nbits, order, false);
            assert_eq!(unpack(&p, nbits, order, false), vals, "nbits={nbits} {order:?}");
        }
    }
}

#[test]
fn msb_and_lsb_differ() {
    let p = [0b0001_1011u8];
    assert_eq!(unpack(&p, 2, BitOrder::MsbFirst, false), vec![0.0, 1.0, 2.0, 3.0]);
    assert_eq!(unpack(&p, 2, BitOrder::LsbFirst, false), vec![3.0, 2.0, 1.0, 0.0]);
}

#[test]
fn sigproc_round_trip() {
    let dir = std::env::temp_dir().join(format!("fbtest_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    for nbits in [2u32, 4, 8, 32] {
        let path = dir.join(format!("t{nbits}.fil"));
        let h = ObsHeader {
            source_name: "J0534+2200".into(),
            telescope: "GMRT".into(),
            backend: "GMRTNEW".into(),
            ra_deg: Some(83.633),
            dec_deg: Some(-0.5),
            az_deg: Some(10.0),
            za_deg: Some(20.0),
            start: Epoch::new(60123, 4567.891),
            tsamp: 81.92e-6,
            nsamples: 0,
            freqs: FreqAxis::uniform(1460.0, -0.5, 16),
            npol: 1,
            pol_type: PolType::Intensity,
            nbits,
            scaling: Scaling::default(),
            barycentric: false,
            beam: 3,
            nbeams: 4,
            order: DataOrder::TimePolFreq,
            provenance: vec![],
        };
        let max = if nbits == 32 { 1000.0 } else { ((1u32 << nbits) - 1) as f32 };
        let vals: Vec<f32> = (0..16 * 40).map(|i| ((i * 7) % 13) as f32 % (max + 1.0)).collect();
        let mut w = SigprocWriter::create(&path, &h, WriteOptions { nbits, data_type: 1, refdm: None, rawdatafile: "x".into() }).unwrap();
        w.write_values(&vals).unwrap();
        w.finish().unwrap();
        let mut r = open(&[&path]).unwrap();
        let rh = r.header();
        assert_eq!(rh.nsamples, 40);
        assert_eq!(rh.nchan(), 16);
        assert!((rh.tsamp - h.tsamp).abs() < 1e-15);
        assert!((rh.start.seconds_since(&h.start)).abs() < 1e-5);
        assert!((rh.dec_deg.unwrap() - -0.5).abs() < 1e-6, "dec sign kept");
        assert_eq!(rh.telescope, "GMRT");
        assert_eq!(rh.beam, 3);
        let b = r.read(0, 40).unwrap();
        assert_eq!(b.data, vals, "nbits={nbits}");
        let b2 = r.read(10, 5).unwrap();
        assert_eq!(&b2.data[..], &vals[160..240]);
    }
}

#[test]
fn dada_reader() {
    use std::io::Write;
    let dir = std::env::temp_dir().join(format!("fbdada_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("x.dada");
    let mut hdr = String::from("HDR_VERSION 1.0\nHDR_SIZE 4096\nNCHAN 8\nNBIT 8\nNPOL 1\nNDIM 1\nTSAMP 64\nBW -16\nFREQ 1400\nUTC_START 2024-01-02-03:04:05\nOBS_OFFSET 0\nSOURCE TESTSRC\nTELESCOPE DAT\n");
    hdr.push_str(&"\0".repeat(4096 - hdr.len()));
    let data: Vec<u8> = (0..8 * 100).map(|i| (i % 251) as u8).collect();
    let mut f = std::fs::File::create(&path).unwrap();
    f.write_all(hdr.as_bytes()).unwrap();
    f.write_all(&data).unwrap();
    drop(f);
    let mut r = open(&[&path]).unwrap();
    let h = r.header();
    assert_eq!(h.nsamples, 100);
    assert_eq!(h.nchan(), 8);
    assert_eq!(h.start.imjd, 60311);
    assert!((h.start.sec - (3.0 * 3600.0 + 4.0 * 60.0 + 5.0)).abs() < 1e-9);
    assert!((h.tsamp - 64e-6).abs() < 1e-15);
    assert!((h.freqs.fch1() - 1407.0).abs() < 1e-9 && h.freqs.is_descending());
    let b = r.read(0, 100).unwrap();
    assert_eq!(b.data, data.iter().map(|v| *v as f32).collect::<Vec<_>>());
}
