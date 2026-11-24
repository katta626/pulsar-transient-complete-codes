use std::f64::consts::PI;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use fb_defs::{Epoch, FbError, Result};
use rustfft::num_complex::Complex64;
use rustfft::FftPlanner;
use serde::{Deserialize, Serialize};

fn dft(x: &[f64]) -> Vec<Complex64> {
    let mut buf: Vec<Complex64> = x.iter().map(|v| Complex64::new(*v, 0.0)).collect();
    FftPlanner::<f64>::new().plan_fft_forward(buf.len()).process(&mut buf);
    buf
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Shift {
    pub phase: f64,
    pub phase_err: f64,
    pub scale: f64,
    pub snr: f64,
}

pub fn fftfit(profile: &[f64], template: &[f64]) -> Result<Shift> {
    let n = profile.len();
    if n != template.len() || n < 8 {
        return Err(FbError::Invalid("profile and template need the same length (>= 8 bins)".into()));
    }
    let p = dft(profile);
    let t = dft(template);
    let nh = n / 2;
    let kmax = nh - 1;
    let a: Vec<f64> = (1..=kmax).map(|k| p[k].norm()).collect();
    let b: Vec<f64> = (1..=kmax).map(|k| t[k].norm()).collect();
    let dphi: Vec<f64> = (1..=kmax).map(|k| t[k].arg() - p[k].arg()).collect();
    let f = |tau: f64| -> f64 { (0..kmax).map(|i| a[i] * b[i] * (dphi[i] - (i + 1) as f64 * tau).cos()).sum() };
    let df = |tau: f64| -> f64 { (0..kmax).map(|i| (i + 1) as f64 * a[i] * b[i] * (dphi[i] - (i + 1) as f64 * tau).sin()).sum() };
    let d2f = |tau: f64| -> f64 { -(0..kmax).map(|i| ((i + 1) as f64).powi(2) * a[i] * b[i] * (dphi[i] - (i + 1) as f64 * tau).cos()).sum::<f64>() };
    let mut best = (0.0, f64::NEG_INFINITY);
    let ngrid = 8 * n;
    for g in 0..ngrid {
        let tau = 2.0 * PI * g as f64 / ngrid as f64;
        let v = f(tau);
        if v > best.1 {
            best = (tau, v);
        }
    }
    let mut tau = best.0;
    for _ in 0..50 {
        let d2 = d2f(tau);
        if d2 >= 0.0 {
            break;
        }
        let step = df(tau) / d2;
        tau -= step;
        if step.abs() < 1e-12 {
            break;
        }
    }
    let sum_b2: f64 = b.iter().map(|v| v * v).sum();
    let scale = f(tau) / sum_b2.max(1e-300);
    let resid: f64 = (0..kmax).map(|i| a[i] * a[i] - 2.0 * scale * a[i] * b[i] * (dphi[i] - (i + 1) as f64 * tau).cos() + scale * scale * b[i] * b[i]).sum();
    let sigma2 = (resid / kmax as f64).max(1e-300);
    let curv = -d2f(tau);
    let tau_err = (sigma2 / (2.0 * scale.abs().max(1e-300) * curv.max(1e-300))).sqrt();
    let phase = (tau / (2.0 * PI)).rem_euclid(1.0);
    let snr = scale * sum_b2.sqrt() / sigma2.sqrt();
    Ok(Shift { phase, phase_err: tau_err / (2.0 * PI), scale, snr })
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Toa {
    pub name: String,
    pub freq_mhz: f64,
    pub epoch: Epoch,
    pub err_us: f64,
    pub site: String,
    pub flags: Vec<(String, String)>,
}

pub fn toa_from_fold(profile: &[f64], template: &[f64], fold_epoch: &Epoch, f0: f64, freq_mhz: f64, site: &str, name: &str) -> Result<Toa> {
    let s = fftfit(profile, template)?;
    let mut ph = s.phase;
    if ph > 0.5 {
        ph -= 1.0;
    }
    let period = 1.0 / f0;
    Ok(Toa {
        name: name.into(),
        freq_mhz,
        epoch: fold_epoch.add_seconds(ph * period),
        err_us: s.phase_err * period * 1e6,
        site: site.into(),
        flags: vec![("snr".into(), format!("{:.2}", s.snr))],
    })
}

pub fn write_tim(path: &Path, toas: &[Toa]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "FORMAT 1")?;
    for t in toas {
        let frac = t.epoch.sec / 86400.0;
        let frac_s = format!("{frac:.16}");
        let flags: String = t.flags.iter().map(|(k, v)| format!(" -{k} {v}")).collect();
        writeln!(w, "{} {:.6} {}.{} {:.3} {}{}", t.name, t.freq_mhz, t.epoch.imjd, &frac_s[2..], t.err_us, t.site, flags)?;
    }
    w.flush()?;
    Ok(())
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Polyco {
    pub psr: String,
    pub tmid: Epoch,
    pub dm: f64,
    pub rphase: f64,
    pub f0: f64,
    pub site: String,
    pub span_min: f64,
    pub obsfreq: f64,
    pub coeffs: Vec<f64>,
}

fn pf(s: &str) -> f64 {
    s.trim().replace(['D', 'd'], "E").parse().unwrap_or(0.0)
}

pub fn read_polycos(path: &Path) -> Result<Vec<Polyco>> {
    let lines: Vec<String> = BufReader::new(File::open(path)?).lines().collect::<std::io::Result<_>>()?;
    let mut out = Vec::new();
    let mut i = 0;
    while i + 1 < lines.len() {
        let l1: Vec<&str> = lines[i].split_whitespace().collect();
        let l2: Vec<&str> = lines[i + 1].split_whitespace().collect();
        if l1.len() < 4 || l2.len() < 6 {
            i += 1;
            continue;
        }
        let tmid_s = l1[3];
        let (ip, fp) = tmid_s.split_once('.').unwrap_or((tmid_s, "0"));
        let tmid = Epoch::new(ip.parse().unwrap_or(0), format!("0.{fp}").parse::<f64>().unwrap_or(0.0) * 86400.0);
        let ncoeff: usize = l2[4].parse().map_err(|_| FbError::Format("bad polyco NCOEFF".into()))?;
        let nlines = ncoeff.div_ceil(3);
        let mut coeffs = Vec::with_capacity(ncoeff);
        for k in 0..nlines {
            let l = lines.get(i + 2 + k).ok_or_else(|| FbError::Format("truncated polyco".into()))?;
            let l = l.trim_end();
            if l.len() >= 25 && l.len() % 25 == 0 {
                for k in 0..l.len() / 25 {
                    coeffs.push(pf(&l[k * 25..(k + 1) * 25]));
                }
            } else {
                let mut tok = String::new();
                let chars: Vec<char> = l.chars().collect();
                for (j, ch) in chars.iter().enumerate() {
                    let prev = if j > 0 { chars[j - 1] } else { ' ' };
                    let starts_new = *ch == '-' && !(prev == 'D' || prev == 'd' || prev == 'E' || prev == 'e' || prev == ' ');
                    if ch.is_whitespace() || starts_new {
                        if !tok.trim().is_empty() {
                            coeffs.push(pf(&tok));
                        }
                        tok.clear();
                        if starts_new {
                            tok.push('-');
                        }
                    } else {
                        tok.push(*ch);
                    }
                }
                if !tok.trim().is_empty() {
                    coeffs.push(pf(&tok));
                }
            }
        }
        coeffs.truncate(ncoeff);
        out.push(Polyco {
            psr: l1[0].into(),
            tmid,
            dm: pf(l1.get(4).copied().unwrap_or("0")),
            rphase: pf(l2[0]),
            f0: pf(l2[1]),
            site: l2[2].into(),
            span_min: pf(l2[3]),
            obsfreq: pf(l2[5]),
            coeffs,
        });
        i += 2 + nlines;
    }
    if out.is_empty() {
        return Err(FbError::Format("no polyco blocks found".into()));
    }
    Ok(out)
}

impl Polyco {
    pub fn covers(&self, t: &Epoch) -> bool {
        (t.seconds_since(&self.tmid) / 60.0).abs() <= self.span_min / 2.0 + 1e-9
    }

    pub fn phase(&self, t: &Epoch) -> f64 {
        let dt = t.seconds_since(&self.tmid) / 60.0;
        let mut poly = 0.0;
        for c in self.coeffs.iter().rev() {
            poly = poly * dt + c;
        }
        self.rphase + dt * 60.0 * self.f0 + poly
    }

    pub fn frequency(&self, t: &Epoch) -> f64 {
        let dt = t.seconds_since(&self.tmid) / 60.0;
        let mut d = 0.0;
        for (i, c) in self.coeffs.iter().enumerate().skip(1).rev() {
            d = d * dt + i as f64 * c;
        }
        self.f0 + d / 60.0
    }
}

pub fn polyco_for<'a>(set: &'a [Polyco], t: &Epoch) -> Option<&'a Polyco> {
    set.iter().filter(|p| p.covers(t)).min_by(|a, b| a.tmid.seconds_since(t).abs().partial_cmp(&b.tmid.seconds_since(t).abs()).unwrap())
}
