use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use fb_defs::{Epoch, ObsHeader, PeriodicCandidate, Result, SinglePulseCandidate};

pub fn write_singlepulse(path: &Path, cands: &[SinglePulseCandidate], h: &ObsHeader) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "# DM      Sigma      Time (s)     Sample    Downfact")?;
    for c in cands {
        let t = c.epoch.seconds_since(&h.start);
        writeln!(w, "{:7.2} {:7.2} {:13.6} {:10} {:3}", c.dm, c.snr, t, c.sample, c.width_samples)?;
    }
    w.flush()?;
    Ok(())
}

pub fn read_singlepulse(path: &Path, h: &ObsHeader) -> Result<Vec<SinglePulseCandidate>> {
    let r = BufReader::new(File::open(path)?);
    let f_ref = h.freqs.mhz.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut out = Vec::new();
    for l in r.lines() {
        let l = l?;
        if l.starts_with('#') || l.trim().is_empty() {
            continue;
        }
        let v: Vec<&str> = l.split_whitespace().collect();
        if v.len() < 5 {
            continue;
        }
        let dm: f64 = v[0].parse().unwrap_or(0.0);
        let snr: f64 = v[1].parse().unwrap_or(0.0);
        let t: f64 = v[2].parse().unwrap_or(0.0);
        let sample: u64 = v[3].parse().unwrap_or(0);
        let w: u32 = v[4].parse().unwrap_or(1);
        out.push(SinglePulseCandidate {
            obs_id: h.observation_id(),
            dm,
            epoch: h.start.add_seconds(t),
            f_ref,
            sample,
            width_samples: w,
            width_s: w as f64 * h.tsamp,
            snr,
            dm_index: 0,
            members: 1,
            beam: h.beam,
            label: None,
            probability: None,
        });
    }
    Ok(out)
}

pub fn read_heimdall(path: &Path, h: &ObsHeader) -> Result<Vec<SinglePulseCandidate>> {
    let r = BufReader::new(File::open(path)?);
    let f_ref = h.freqs.mhz.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut out = Vec::new();
    for l in r.lines() {
        let l = l?;
        let v: Vec<&str> = l.split_whitespace().collect();
        if v.len() < 7 {
            continue;
        }
        let snr: f64 = v[0].parse().unwrap_or(0.0);
        let sample: u64 = v[1].parse().unwrap_or(0);
        let t: f64 = v[2].parse().unwrap_or(sample as f64 * h.tsamp);
        let filt: u32 = v[3].parse().unwrap_or(0);
        let dm_index: u32 = v[4].parse().unwrap_or(0);
        let dm: f64 = v[5].parse().unwrap_or(0.0);
        let members: u32 = v[6].parse().unwrap_or(1);
        let width = 1u32 << filt;
        out.push(SinglePulseCandidate {
            obs_id: h.observation_id(),
            dm,
            epoch: h.start.add_seconds(t),
            f_ref,
            sample,
            width_samples: width,
            width_s: width as f64 * h.tsamp,
            snr,
            dm_index,
            members,
            beam: h.beam,
            label: None,
            probability: None,
        });
    }
    Ok(out)
}

pub fn write_heimdall(path: &Path, cands: &[SinglePulseCandidate], h: &ObsHeader) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    for c in cands {
        let filt = (c.width_samples.max(1) as f64).log2().round() as u32;
        let t = c.epoch.seconds_since(&h.start);
        writeln!(w, "{:.6}\t{}\t{:.6}\t{}\t{}\t{:.4}\t{}\t{}\t{}", c.snr, c.sample, t, filt, c.dm_index, c.dm, c.members, c.sample, c.sample + c.width_samples as u64)?;
    }
    w.flush()?;
    Ok(())
}

pub fn write_csv_single(path: &Path, cands: &[SinglePulseCandidate]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "obs_id,dm,mjd,f_ref_mhz,sample,width_samples,width_s,snr,dm_index,members,beam,label,probability")?;
    for c in cands {
        writeln!(
            w,
            "{},{:.6},{:.12},{:.6},{},{},{:.9},{:.4},{},{},{},{},{}",
            c.obs_id,
            c.dm,
            c.epoch.mjd(),
            c.f_ref,
            c.sample,
            c.width_samples,
            c.width_s,
            c.snr,
            c.dm_index,
            c.members,
            c.beam,
            c.label.clone().unwrap_or_default(),
            c.probability.map(|p| format!("{p:.4}")).unwrap_or_default()
        )?;
    }
    w.flush()?;
    Ok(())
}

pub fn write_csv_periodic(path: &Path, cands: &[PeriodicCandidate]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "obs_id,f_hz,period_s,fdot,fddot,pdot,epoch_mjd,dm,snr,sigma,numharm,method,label,probability")?;
    for c in cands {
        writeln!(
            w,
            "{},{:.12},{:.12},{:.6e},{:.6e},{:.6e},{:.12},{:.4},{:.3},{:.3},{},{},{},{}",
            c.obs_id,
            c.f_hz,
            c.period(),
            c.fdot,
            c.fddot,
            c.pdot(),
            c.epoch.mjd(),
            c.dm,
            c.snr,
            c.sigma,
            c.numharm,
            c.method,
            c.label.clone().unwrap_or_default(),
            c.probability.map(|p| format!("{p:.4}")).unwrap_or_default()
        )?;
    }
    w.flush()?;
    Ok(())
}

pub fn epoch_from_mjd(mjd: f64) -> Epoch {
    Epoch::from_mjd(mjd)
}
