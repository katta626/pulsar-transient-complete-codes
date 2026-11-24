use std::path::{Path, PathBuf};

use fb_defs::{Epoch, FbError, Result};
use fb_timing::{toa_from_fold, write_tim};
use serde_json::Value;

use crate::args::Args;

pub const TOA_USAGE: &str = "fb toa FOLD.json [--template FOLD_OR_PROFILE.json] [--freq MHz] [--site CODE] [-o OUT.tim]";

fn load(p: &Path) -> Result<Value> {
    Ok(serde_json::from_str(&std::fs::read_to_string(p)?)?)
}

fn arr(v: &Value, k: &str) -> Result<Vec<f64>> {
    v.get(k)
        .and_then(|x| x.as_array())
        .map(|a| a.iter().filter_map(|y| y.as_f64()).collect())
        .ok_or_else(|| FbError::Format(format!("fold file has no '{k}'")))
}

pub fn toa(raw: &[String]) -> Result<()> {
    let a = Args::parse(raw, &["template", "freq", "site", "o"]).map_err(FbError::Invalid)?;
    let p = PathBuf::from(a.positional.first().ok_or_else(|| FbError::Invalid(format!("usage: {TOA_USAGE}")))?);
    let f = load(&p)?;
    let nsub = f["nsub"].as_u64().unwrap_or(1) as usize;
    let nchan = f["nchan"].as_u64().unwrap_or(1) as usize;
    let nbins = f["nbins"].as_u64().ok_or_else(|| FbError::Format("no nbins".into()))? as usize;
    let data = arr(&f, "data")?;
    let counts = arr(&f, "counts")?;
    let period = f["period"].as_f64().ok_or_else(|| FbError::Format("no period".into()))?;
    let pdot = f["pdot"].as_f64().unwrap_or(0.0);
    let epoch0 = Epoch::from_mjd(f["epoch_mjd"].as_f64().ok_or_else(|| FbError::Format("no epoch_mjd".into()))?);
    let template = match a.get("template") {
        Some(t) => arr(&load(Path::new(t))?, "profile")?,
        None => arr(&f, "profile")?,
    };
    if template.len() != nbins {
        return Err(FbError::Incompatible(format!("template has {} bins, fold has {nbins}", template.len())));
    }
    let freq: f64 = a.num("freq", 1400.0).map_err(FbError::Invalid)?;
    let site = a.get("site").unwrap_or("@").to_string();
    let f0 = 1.0 / period;
    let fdot = -pdot / (period * period);
    let mut toas = Vec::new();
    for s in 0..nsub {
        let mut prof = vec![0f64; nbins];
        let mut cnt = vec![0f64; nbins];
        for c in 0..nchan {
            for b in 0..nbins {
                let k = (s * nchan + c) * nbins + b;
                prof[b] += data[k];
                cnt[b] += counts[k];
            }
        }
        if cnt.iter().all(|x| *x == 0.0) {
            continue;
        }
        let prof: Vec<f64> = prof.iter().zip(&cnt).map(|(v, n)| if *n > 0.0 { v / n } else { 0.0 }).collect();
        let tsub = f["duration_s"].as_f64().unwrap_or(0.0) * (s as f64 + 0.5) / nsub as f64;
        let n_turns = (f0 * tsub + 0.5 * fdot * tsub * tsub).floor();
        let t_k = if fdot == 0.0 { n_turns / f0 } else { tsub - ((f0 * tsub + 0.5 * fdot * tsub * tsub) - n_turns) / (f0 + fdot * tsub) };
        let ep = epoch0.add_seconds(t_k);
        let name = format!("{}_{s:04}", p.file_stem().and_then(|x| x.to_str()).unwrap_or("fold"));
        toas.push(toa_from_fold(&prof, &template, &ep, f0 + fdot * t_k, freq, &site, &name)?);
    }
    let out = a.get("o").map(PathBuf::from).unwrap_or_else(|| p.with_extension("tim"));
    write_tim(&out, &toas)?;
    println!("{} TOAs -> {}", toas.len(), out.display());
    for t in toas.iter().take(5) {
        println!("  {} {:.12} ± {:.2} us", t.name, t.epoch.mjd(), t.err_us);
    }
    Ok(())
}
