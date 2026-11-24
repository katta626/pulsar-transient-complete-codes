use std::path::{Path, PathBuf};

use fb_defs::{FbError, Result};
use fb_io::{open_with, PsrfitsOptions, ScalePolicy, SigprocWriter, WriteOptions, ZeroOffPolicy};

use crate::args::Args;

pub const USAGE: &str = "fb convert INPUT... [-o OUT.fil] [--nbits N] [--scales auto|always|never]
           [--zero-off spec|dspsr|ignore] [--keep-pols] [--order descending|keep]
           [--gulp SAMPLES] [--max-gap SECONDS] [--overwrite]";

pub fn run(raw: &[String]) -> Result<()> {
    let a = Args::parse(raw, &["o", "output", "nbits", "scales", "zero-off", "order", "gulp", "max-gap"]).map_err(FbError::Invalid)?;
    if a.positional.is_empty() {
        return Err(FbError::Invalid(format!("usage: {USAGE}")));
    }
    let inputs: Vec<PathBuf> = a.positional.iter().map(PathBuf::from).collect();
    let out = a
        .get("o")
        .or(a.get("output"))
        .map(PathBuf::from)
        .unwrap_or_else(|| inputs[0].with_extension("fil"));
    if out.exists() && !a.flag("overwrite") {
        return Err(FbError::Invalid(format!("{} exists; pass --overwrite", out.display())));
    }

    let scales = match a.get("scales").unwrap_or("auto") {
        "auto" => ScalePolicy::Auto,
        "always" => ScalePolicy::Always,
        "never" => ScalePolicy::Never,
        s => return Err(FbError::Invalid(format!("--scales {s}"))),
    };
    let zero_off_req = a.get("zero-off");
    let sum_pols = !a.flag("keep-pols");

    let probe = open_with(&inputs, PsrfitsOptions { zero_off: ZeroOffPolicy::Ignore, scales, weights: true, sum_pols: false })?;
    let ph = probe.header();
    let scaled = ph.provenance.iter().any(|p| p.contains("scales=true"));
    drop(probe);

    let lossless = !scaled && zero_off_req.is_none();
    let zero_off = match zero_off_req {
        Some("spec") => ZeroOffPolicy::Spec,
        Some("dspsr") => ZeroOffPolicy::Dspsr,
        Some("ignore") | None if lossless => ZeroOffPolicy::Ignore,
        Some("ignore") => ZeroOffPolicy::Ignore,
        None => ZeroOffPolicy::Spec,
        Some(s) => return Err(FbError::Invalid(format!("--zero-off {s}"))),
    };
    let opts = PsrfitsOptions { zero_off, scales, weights: true, sum_pols };
    let mut reader = open_with(&inputs, opts)?;
    let h = reader.header();
    let max_gap: f64 = a.num("max-gap", 60.0f64).map_err(FbError::Invalid)?;
    let total_gap: f64 = reader.gaps().iter().map(|g| g.len as f64 * h.tsamp).sum();
    if total_gap > max_gap {
        return Err(FbError::Invalid(format!(
            "inputs contain {:.1} s of missing data (limit {max_gap} s); convert the files separately or raise --max-gap",
            total_gap
        )));
    }

    let widen = sum_pols && ph.npol > 1;
    let default_bits = if lossless && !widen { h.nbits } else { 32 };
    let nbits: u32 = a.num("nbits", default_bits).map_err(FbError::Invalid)?;
    let descending = a.get("order").unwrap_or("descending") == "descending";
    let flip = descending && !h.freqs.is_descending();

    let mut out_h = (*h).clone();
    if flip {
        out_h.freqs = h.freqs.reversed();
    }
    out_h.nbits = nbits;
    out_h.provenance.push(format!("convert nbits={nbits} flip={flip} lossless={lossless}"));
    let raw_name = inputs[0].file_name().and_then(|s| s.to_str()).unwrap_or("").to_string();
    let mut w = SigprocWriter::create(&out, &out_h, WriteOptions { nbits, data_type: 1, refdm: None, rawdatafile: raw_name })?;

    let gulp: usize = a.num("gulp", 65536usize).map_err(FbError::Invalid)?;
    let nchan = h.nchan();
    let mut s = 0u64;
    while s < h.nsamples {
        let mut b = reader.read(s, gulp)?;
        if flip {
            for spec in b.data.chunks_mut(nchan) {
                spec.reverse();
            }
        }
        w.write_block(&b)?;
        s += b.nsamp as u64;
        if b.nsamp == 0 {
            break;
        }
    }
    let n = w.finish()?;
    for g in reader.gaps() {
        eprintln!("gap filled: samples {}..{} ({} s)", g.start, g.start + g.len, g.len as f64 * h.tsamp);
    }
    println!("Wrote {} ({} samples, {} ch, {}-bit{})", out.display(), n, nchan, nbits, if flip { ", band flipped to descending" } else { "" });
    Ok(())
}

pub fn header(path: &Path) -> Result<()> {
    let r = open_with(&[path], PsrfitsOptions::raw())?;
    let h = r.header();
    println!("{}", serde_json::to_string_pretty(&*h).map_err(|e| FbError::Format(e.to_string()))?);
    for g in r.gaps() {
        println!("gap: start={} len={}", g.start, g.len);
    }
    Ok(())
}
