use std::path::{Path, PathBuf};

use fb_defs::{FbError, Result, RfiMask};
use fb_formats::candidates::{write_csv_single, write_heimdall, write_singlepulse};
use fb_formats::presto::PrestoMask;
use fb_formats::store::{load_mask, save_json, save_mask};
use fb_io::{open_with, PsrfitsOptions};
use fb_pipeline::{run_single_pulse, Dedisperser, SinglePulseConfig};

use crate::args::Args;

pub const SP_USAGE: &str = "fb sp INPUT... [-o PREFIX] [--dm LO:HI] [--dm-step S] [--dm-tol 1.25] [--width-s 0.001]
       [--threshold 6] [--max-width 4096] [--gulp N] [--subband NSUB] [--mask FILE]
       [--zero-dm] [--no-iqrm] [--sk] [--stats] [--format csv,presto,heimdall]";

pub const RFI_USAGE: &str = "fb rfi INPUT... [-o MASK.json] [--interval 1.0] [--iqrm-threshold 3] [--sk] [--stats]
       [--presto OUT.mask]";

fn load_any_mask(p: &Path, h: &fb_defs::ObsHeader) -> Result<RfiMask> {
    if p.extension().and_then(|e| e.to_str()) == Some("mask") {
        Ok(PrestoMask::read(p)?.to_mask_for(h))
    } else {
        load_mask(p)
    }
}

fn config(a: &Args) -> Result<SinglePulseConfig> {
    let e = |s: String| FbError::Invalid(s);
    let mut cfg = SinglePulseConfig::default();
    if let Some(r) = a.get("dm") {
        let (lo, hi) = r.split_once(':').ok_or_else(|| FbError::Invalid("--dm LO:HI".into()))?;
        cfg.dm_lo = lo.parse().map_err(|_| FbError::Invalid("--dm".into()))?;
        cfg.dm_hi = hi.parse().map_err(|_| FbError::Invalid("--dm".into()))?;
    }
    if a.get("dm-step").is_some() {
        cfg.dm_step = Some(a.num("dm-step", 1.0).map_err(e)?);
    }
    cfg.dm_tol = a.num("dm-tol", 1.25).map_err(e)?;
    cfg.pulse_width = a.num("width-s", 1e-3).map_err(e)?;
    cfg.search.threshold = a.num("threshold", 6.0f32).map_err(e)?;
    cfg.search.max_width = a.num("max-width", 4096usize).map_err(e)?;
    cfg.gulp = a.num("gulp", 1usize << 16).map_err(e)?;
    if a.get("subband").is_some() {
        cfg.dedisperser = Dedisperser::Subband(a.num("subband", 32usize).map_err(e)?);
    }
    cfg.rfi.zero_dm = a.flag("zero-dm");
    cfg.rfi.iqrm = !a.flag("no-iqrm");
    cfg.rfi.sk = a.flag("sk");
    cfg.rfi.interval_stats = a.flag("stats");
    cfg.rfi.interval = a.num("interval", 1.0).map_err(e)?;
    cfg.rfi.iqrm_threshold = a.num("iqrm-threshold", 3.0).map_err(e)?;
    Ok(cfg)
}

const VALUED: &[&str] = &["o", "dm", "dm-step", "dm-tol", "width-s", "threshold", "max-width", "gulp", "subband", "mask", "format", "interval", "iqrm-threshold", "presto"];

pub fn sp(raw: &[String]) -> Result<()> {
    let a = Args::parse(raw, VALUED).map_err(FbError::Invalid)?;
    if a.positional.is_empty() {
        return Err(FbError::Invalid(format!("usage: {SP_USAGE}")));
    }
    let inputs: Vec<PathBuf> = a.positional.iter().map(PathBuf::from).collect();
    let mut reader = open_with(&inputs, PsrfitsOptions::default())?;
    let h = reader.header();
    let mut cfg = config(&a)?;
    if let Some(m) = a.get("mask") {
        cfg.external_mask = Some(load_any_mask(Path::new(m), &h)?);
    }
    let prefix = a.get("o").map(PathBuf::from).unwrap_or_else(|| inputs[0].with_extension(""));
    let res = run_single_pulse(reader.as_mut(), &cfg)?;
    let formats = a.get("format").unwrap_or("csv,presto");
    for f in formats.split(',') {
        match f.trim() {
            "csv" => write_csv_single(&prefix.with_extension("cands.csv"), &res.candidates)?,
            "presto" => write_singlepulse(&prefix.with_extension("singlepulse"), &res.candidates, &h)?,
            "heimdall" => write_heimdall(&prefix.with_extension("cand"), &res.candidates, &h)?,
            "" => {}
            other => return Err(FbError::Invalid(format!("unknown format '{other}'"))),
        }
    }
    save_mask(&prefix.with_extension("mask.json"), &res.mask)?;
    save_json(&prefix.with_extension("dmplan.json"), &res.plan)?;
    println!(
        "{} candidates, {} DM trials ({:.1}..{:.1}), {:.1}% of the mask flagged -> {}.*",
        res.candidates.len(),
        res.plan.ndm(),
        cfg.dm_lo,
        cfg.dm_hi,
        100.0 * res.mask.fraction_bad(),
        prefix.display()
    );
    Ok(())
}

pub fn rfi(raw: &[String]) -> Result<()> {
    let a = Args::parse(raw, VALUED).map_err(FbError::Invalid)?;
    if a.positional.is_empty() {
        return Err(FbError::Invalid(format!("usage: {RFI_USAGE}")));
    }
    let inputs: Vec<PathBuf> = a.positional.iter().map(PathBuf::from).collect();
    let mut reader = open_with(&inputs, PsrfitsOptions::default())?;
    let h = reader.header();
    let cfg = config(&a)?;
    let mut mask = RfiMask::empty(&h, cfg.rfi.interval, "filterbank");
    let per = ((cfg.rfi.interval / h.tsamp).round() as usize).max(1);
    let gulp = per * ((65536 / per).max(1));
    let mut s = 0u64;
    while s < h.nsamples {
        let mut b = reader.read(s, gulp)?;
        if b.nsamp == 0 {
            break;
        }
        if cfg.rfi.zero_dm {
            fb_rfi::zero_dm(&mut b, true)?;
        }
        if cfg.rfi.interval_stats {
            fb_rfi::interval_stats_mask(&b, cfg.rfi.interval_threshold, &mut mask)?;
        }
        if cfg.rfi.sk {
            fb_rfi::spectral_kurtosis(&b, cfg.rfi.sk_m, 1.0, cfg.rfi.sk_sigma, &mut mask)?;
        }
        if cfg.rfi.iqrm {
            let r = (b.nchan() / 10).max(1);
            fb_rfi::iqrm_mask(&b, fb_rfi::Statistic::Std, r, cfg.rfi.iqrm_threshold, &mut mask)?;
        }
        s += b.nsamp as u64;
    }
    let out = a.get("o").map(PathBuf::from).unwrap_or_else(|| inputs[0].with_extension("mask.json"));
    save_mask(&out, &mask)?;
    if let Some(p) = a.get("presto") {
        PrestoMask::from_mask(&mask, h.tsamp)?.write(Path::new(p))?;
    }
    println!("{:.1}% flagged -> {}", 100.0 * mask.fraction_bad(), out.display());
    Ok(())
}

pub fn export_mask(raw: &[String]) -> Result<()> {
    let a = Args::parse(raw, VALUED).map_err(FbError::Invalid)?;
    if a.positional.len() != 2 {
        return Err(FbError::Invalid("usage: fb export-mask DATAFILE MASK(.json|.mask) -o OUT(.json|.mask)".into()));
    }
    let r = open_with(&[&a.positional[0]], PsrfitsOptions::default())?;
    let h = r.header();
    let m = load_any_mask(Path::new(&a.positional[1]), &h)?;
    let out = PathBuf::from(a.get("o").ok_or_else(|| FbError::Invalid("-o OUT".into()))?);
    if out.extension().and_then(|e| e.to_str()) == Some("mask") {
        PrestoMask::from_mask(&m, h.tsamp)?.write(&out)?;
    } else {
        save_mask(&out, &m)?;
    }
    println!("{:.1}% flagged -> {}", 100.0 * m.fraction_bad(), out.display());
    Ok(())
}
