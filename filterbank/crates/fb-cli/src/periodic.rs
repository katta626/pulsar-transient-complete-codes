use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use fb_defs::{Epoch, FbError, PeriodicCandidate, Result, Series};
use fb_formats::candidates::write_csv_periodic;
use fb_formats::presto::{read_dat, PrestoMask};
use fb_formats::store::{load_mask, save_json};
use fb_io::{open_with, PsrfitsOptions};
use fb_pipeline::{run_periodic, PeriodicConfig};
use serde_json::json;

use crate::args::Args;

pub const PERIODIC_USAGE: &str = "fb periodic INPUT... [-o PREFIX] [--dm LO:HI] [--dm-step S] [--fmin 0.1] [--fmax 2000]
       [--harm 16] [--sigma 6] [--accel MAX] [--accel-tol 1.1] [--zmax Z] [--wmax W] [--dz 2] [--dw 20] [--refine 50] [--ffa] [--pmin 0.5] [--pmax 30]
       [--no-fft] [--fold N] [--nbins 64] [--dedisp auto|brute|subband|fdmt] [--min-dm-hits 1] [--zap FREQ:WIDTH,...] [--mask FILE] [--zero-dm]";

pub const SERIES_USAGE: &str = "fb fft SERIES(.dat|.tim) [--fmin] [--fmax] [--harm] [--sigma] [--accel MAX | --zmax Z [--wmax W]] [-o OUT.csv]
fb ffa SERIES(.dat|.tim) [--pmin] [--pmax] [--snr 7] [-o OUT.csv]
fb fold SERIES|FILTERBANK --period P [--pdot PD] [--dm DM] [--nbins 64] [-o OUT.json]";

const VALUED: &[&str] = &["o", "dm", "dm-step", "fmin", "fmax", "harm", "sigma", "accel", "accel-tol", "pmin", "pmax", "fold", "nbins", "min-dm-hits", "zap", "mask", "snr", "period", "pdot", "nsub", "gulp", "dedisp", "subband", "zmax", "wmax", "dz", "dw", "refine"];

fn e(s: String) -> FbError {
    FbError::Invalid(s)
}

fn fft_params(a: &Args) -> Result<fb_fft::FftSearchParams> {
    Ok(fb_fft::FftSearchParams {
        f_min: a.num("fmin", 0.1).map_err(e)?,
        f_max: a.num("fmax", 2000.0).map_err(e)?,
        max_harm: a.num("harm", 16u32).map_err(e)?,
        sigma_threshold: a.num("sigma", 6.0).map_err(e)?,
        max_cands: 500,
    })
}

fn fdas_params(a: &Args) -> Result<Option<fb_fft::fdas::FdasParams>> {
    if a.get("zmax").is_none() && a.get("wmax").is_none() {
        return Ok(None);
    }
    Ok(Some(fb_fft::fdas::FdasParams {
        zmax: a.num("zmax", 0.0).map_err(e)?,
        dz: a.num("dz", 2.0).map_err(e)?,
        wmax: a.num("wmax", 0.0).map_err(e)?,
        dw: a.num("dw", 20.0).map_err(e)?,
        refine: a.num("refine", 50usize).map_err(e)?,
    }))
}

fn ffa_params(a: &Args) -> Result<fb_ffa::FfaParams> {
    Ok(fb_ffa::FfaParams {
        period_min: a.num("pmin", 0.5).map_err(e)?,
        period_max: a.num("pmax", 30.0).map_err(e)?,
        snr_threshold: a.num("snr", 7.0).map_err(e)?,
        ..Default::default()
    })
}

fn birdies(a: &Args) -> Result<Vec<(f64, f64)>> {
    let mut out = Vec::new();
    if let Some(z) = a.get("zap") {
        for item in z.split(',') {
            let (f, w) = item.split_once(':').unwrap_or((item, "0.05"));
            out.push((f.parse().map_err(|_| e(format!("--zap {item}")))?, w.parse().map_err(|_| e(format!("--zap {item}")))?));
        }
    }
    Ok(out)
}

pub fn periodic(raw: &[String]) -> Result<()> {
    let a = Args::parse(raw, VALUED).map_err(e)?;
    if a.positional.is_empty() {
        return Err(e(format!("usage: {PERIODIC_USAGE}")));
    }
    let inputs: Vec<PathBuf> = a.positional.iter().map(PathBuf::from).collect();
    let mut reader = open_with(&inputs, PsrfitsOptions::default())?;
    let h = reader.header();
    let mut cfg = PeriodicConfig::default();
    if let Some(r) = a.get("dm") {
        let (lo, hi) = r.split_once(':').ok_or_else(|| e("--dm LO:HI".into()))?;
        cfg.dm_lo = lo.parse().map_err(|_| e("--dm".into()))?;
        cfg.dm_hi = hi.parse().map_err(|_| e("--dm".into()))?;
    }
    if a.get("dm-step").is_some() {
        cfg.dm_step = Some(a.num("dm-step", 1.0).map_err(e)?);
    }
    cfg.fft = if a.flag("no-fft") { None } else { Some(fft_params(&a)?) };
    cfg.max_accel = a.num("accel", 0.0).map_err(e)?;
    cfg.accel_tol = a.num("accel-tol", 1.1).map_err(e)?;
    cfg.fdas = fdas_params(&a)?;
    cfg.ffa = if a.flag("ffa") { Some(ffa_params(&a)?) } else { None };
    cfg.fold_top = a.num("fold", 20usize).map_err(e)?;
    cfg.nbins = a.num("nbins", 64usize).map_err(e)?;
    cfg.min_dm_hits = a.num("min-dm-hits", 1usize).map_err(e)?;
    cfg.birdies = birdies(&a)?;
    cfg.rfi.zero_dm = a.flag("zero-dm");
    cfg.dedisperser = fb_pipeline::Dedisperser::parse(a.get("dedisp").unwrap_or("auto"), a.num("subband", 64usize).map_err(e)?)
        .ok_or_else(|| e("--dedisp auto|brute|subband|fdmt".into()))?;
    cfg.gulp = a.num("gulp", 1usize << 18).map_err(e)?;
    if let Some(m) = a.get("mask") {
        let p = Path::new(m);
        cfg.external_mask = Some(if p.extension().and_then(|x| x.to_str()) == Some("mask") { PrestoMask::read(p)?.to_mask_for(&h) } else { load_mask(p)? });
    }
    let prefix = a.get("o").map(PathBuf::from).unwrap_or_else(|| inputs[0].with_extension(""));
    let res = run_periodic(reader.as_mut(), &cfg)?;
    write_csv_periodic(&fb_defs::add_ext(&prefix, "periodic.csv"), &res.sifted)?;
    let folds: Vec<_> = res
        .folded
        .iter()
        .map(|f| json!({"candidate": f.candidate, "reduced_chi2": f.reduced_chi2, "profile": f.profile, "nsub": f.folded.nsub, "nchan": f.folded.nchan, "nbins": f.folded.nbins}))
        .collect();
    save_json(&fb_defs::add_ext(&prefix, "folds.json"), &folds)?;
    println!("{} raw -> {} sifted candidates over {} DMs; folded {}", res.raw.len(), res.sifted.len(), res.plan.ndm(), res.folded.len());
    for f in res.folded.iter().take(10) {
        println!(
            "  P={:.6} s  DM={:.2}  sigma={:.1}  chi2={:.2}  {}",
            f.candidate.period(),
            f.candidate.dm,
            f.candidate.sigma,
            f.reduced_chi2,
            f.candidate.method
        );
    }
    Ok(())
}

fn load_series(p: &Path) -> Result<Series> {
    if p.extension().and_then(|x| x.to_str()) == Some("dat") || fb_defs::add_ext(p, "inf").exists() && p.extension().and_then(|x| x.to_str()) != Some("fil") {
        return read_dat(p);
    }
    let mut r = open_with(&[p], PsrfitsOptions::default())?;
    let h = r.header();
    if h.nchan() != 1 {
        return Err(e("input has many channels; dedisperse first or use `fb periodic`".into()));
    }
    let b = r.read(0, h.nsamples as usize)?;
    let f_ref = h.freqs.fch1();
    Ok(Series { header: (*h).clone(), dm: 0.0, f_ref, data: b.data })
}

fn print_top(c: &[PeriodicCandidate]) {
    for x in c.iter().take(10) {
        println!("  P={:.9} s  f={:.6} Hz  sigma={:.2}  snr={:.2}  {}", x.period(), x.f_hz, x.sigma, x.snr, x.method);
    }
}

pub fn fft(raw: &[String]) -> Result<()> {
    let a = Args::parse(raw, VALUED).map_err(e)?;
    let p = a.positional.first().ok_or_else(|| e(format!("usage: {SERIES_USAGE}")))?;
    let s = load_series(Path::new(p))?;
    let t_obs = s.len() as f64 * s.tsamp();
    let trials = fb_fft::accel_trials(t_obs, s.tsamp(), a.num("accel", 0.0).map_err(e)?, 1e-3, a.num("accel-tol", 1.1).map_err(e)?);
    let c = match fdas_params(&a)? {
        Some(fp) => fb_fft::fdas::search_series(&s, &birdies(&a)?, &fp, &fft_params(&a)?),
        None => fb_fft::search_series(&s, &trials, &birdies(&a)?, &fft_params(&a)?),
    };
    let sifted = fb_fft::sift(&c, t_obs, 1, 16);
    let out = a.get("o").map(PathBuf::from).unwrap_or_else(|| Path::new(p).with_extension("fft.csv"));
    write_csv_periodic(&out, &sifted)?;
    println!("{} candidates ({} accel trials) -> {}", sifted.len(), trials.len(), out.display());
    print_top(&sifted);
    Ok(())
}

pub fn ffa(raw: &[String]) -> Result<()> {
    let a = Args::parse(raw, VALUED).map_err(e)?;
    let p = a.positional.first().ok_or_else(|| e(format!("usage: {SERIES_USAGE}")))?;
    let s = load_series(Path::new(p))?;
    let c = fb_ffa::search(&s, &ffa_params(&a)?);
    let sifted = fb_fft::sift(&c, s.len() as f64 * s.tsamp(), 1, 16);
    let out = a.get("o").map(PathBuf::from).unwrap_or_else(|| Path::new(p).with_extension("ffa.csv"));
    write_csv_periodic(&out, &sifted)?;
    println!("{} candidates -> {}", sifted.len(), out.display());
    print_top(&sifted);
    Ok(())
}

pub fn fold(raw: &[String]) -> Result<()> {
    let a = Args::parse(raw, VALUED).map_err(e)?;
    let p = PathBuf::from(a.positional.first().ok_or_else(|| e(format!("usage: {SERIES_USAGE}")))?);
    let period: f64 = a.get("period").ok_or_else(|| e("--period P".into()))?.parse().map_err(|_| e("--period".into()))?;
    let pdot: f64 = a.num("pdot", 0.0).map_err(e)?;
    let dm: f64 = a.num("dm", 0.0).map_err(e)?;
    let nbins: usize = a.num("nbins", 64usize).map_err(e)?;
    let nsub: usize = a.num("nsub", 32usize).map_err(e)?;
    let f = 1.0 / period;
    let fdot = -pdot / (period * period);
    let is_series = p.extension().and_then(|x| x.to_str()) == Some("dat") || p.extension().and_then(|x| x.to_str()) == Some("tim");
    let duration;
    let (epoch, folded, chi): (Epoch, fb_fold::Folded, f64) = if is_series {
        let s = load_series(&p)?;
        duration = s.len() as f64 * s.tsamp();
        let fp = fb_fold::FoldParams { f, fdot, fddot: 0.0, epoch: s.header.start, dm, nbins };
        let (_, sd) = s.mean_std();
        let fo = fb_fold::fold_series(&s, &fp, nsub);
        let chi = fo.reduced_chi2(sd * sd);
        (s.header.start, fo, chi)
    } else {
        let mut r = open_with(&[&p], PsrfitsOptions::default())?;
        let h = r.header();
        duration = h.duration();
        let c = PeriodicCandidate {
            obs_id: h.observation_id(),
            f_hz: f,
            fdot,
            fddot: 0.0,
            epoch: h.start,
            dm,
            snr: 0.0,
            sigma: 0.0,
            numharm: 0,
            method: "manual".into(),
            label: None,
            probability: None,
        };
        let mask = fb_defs::RfiMask::empty(&h, 1.0, "none");
        let fc = fb_pipeline::periodic::fold_candidate(r.as_mut(), &c, &mask, nbins, nsub, 16, 1 << 18)?;
        (h.start, fc.folded, fc.reduced_chi2)
    };
    let out = a.get("o").map(PathBuf::from).unwrap_or_else(|| p.with_extension("fold.json"));
    let prof = fb_fold::normalised_profile(&folded);
    let mut w = BufWriter::new(File::create(&out)?);
    let doc = json!({"period": period, "pdot": pdot, "dm": dm, "epoch_mjd": epoch.mjd(), "duration_s": duration, "reduced_chi2": chi, "nbins": nbins,
        "nsub": folded.nsub, "nchan": folded.nchan, "profile": prof, "data": folded.data, "counts": folded.counts});
    serde_json::to_writer(&mut w, &doc).map_err(|x| FbError::Format(x.to_string()))?;
    w.flush()?;
    println!("reduced chi2 = {chi:.2} -> {}", out.display());
    Ok(())
}
