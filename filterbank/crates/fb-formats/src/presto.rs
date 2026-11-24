use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use fb_defs::{DataOrder, Epoch, FbError, FreqAxis, ObsHeader, PolType, Result, RfiMask, Scaling, Series};

#[derive(Debug, Clone, Default)]
pub struct InfData {
    pub basename: String,
    pub telescope: String,
    pub instrument: String,
    pub object: String,
    pub ra: String,
    pub dec: String,
    pub observer: String,
    pub epoch: Epoch,
    pub barycentric: bool,
    pub n: u64,
    pub dt: f64,
    pub breaks: Vec<(u64, u64)>,
    pub band: String,
    pub beam_diam: f64,
    pub dm: f64,
    pub lofreq: f64,
    pub bw: f64,
    pub numchan: usize,
    pub chan_wid: f64,
    pub analyzer: String,
    pub notes: String,
}

fn line(w: &mut impl Write, label: &str, value: impl std::fmt::Display) -> std::io::Result<()> {
    writeln!(w, " {:<38} =  {}", label, value)
}

pub fn root_of(path: &Path) -> PathBuf {
    let mut p = path.to_path_buf();
    if matches!(p.extension().and_then(|e| e.to_str()), Some("dat" | "inf" | "fft" | "mask")) {
        p.set_extension("");
    }
    p
}

impl InfData {
    pub fn write(&self, path: &Path) -> Result<()> {
        let mut w = BufWriter::new(File::create(path)?);
        line(&mut w, "Data file name without suffix", &self.basename)?;
        line(&mut w, "Telescope used", &self.telescope)?;
        line(&mut w, "Instrument used", &self.instrument)?;
        line(&mut w, "Object being observed", &self.object)?;
        line(&mut w, "J2000 Right Ascension (hh:mm:ss.ssss)", &self.ra)?;
        line(&mut w, "J2000 Declination     (dd:mm:ss.ssss)", &self.dec)?;
        line(&mut w, "Data observed by", &self.observer)?;
        let frac = format!("{:.15}", self.epoch.sec / 86400.0);
        line(&mut w, "Epoch of observation (MJD)", format!("{}.{}", self.epoch.imjd, frac.trim_start_matches("0.")))?;
        line(&mut w, "Barycentered?           (1 yes, 0 no)", self.barycentric as i32)?;
        line(&mut w, "Number of bins in the time series", self.n)?;
        line(&mut w, "Width of each time series bin (sec)", format!("{:.15e}", self.dt))?;
        line(&mut w, "Any breaks in the data? (1 yes, 0 no)", (!self.breaks.is_empty()) as i32)?;
        for (i, (a, b)) in self.breaks.iter().enumerate() {
            writeln!(w, " On/Off bin pair #{:>3}                 =  {:<11}, {:<11}", i + 1, a, b)?;
        }
        line(&mut w, "Type of observation (EM band)", if self.band.is_empty() { "Radio" } else { &self.band })?;
        line(&mut w, "Beam diameter (arcsec)", format!("{:.0}", self.beam_diam))?;
        line(&mut w, "Dispersion measure (cm-3 pc)", format!("{:.12}", self.dm))?;
        line(&mut w, "Central freq of low channel (MHz)", format!("{:.12}", self.lofreq))?;
        line(&mut w, "Total bandwidth (MHz)", format!("{:.12}", self.bw))?;
        line(&mut w, "Number of channels", self.numchan)?;
        line(&mut w, "Channel bandwidth (MHz)", format!("{:.12}", self.chan_wid))?;
        line(&mut w, "Data analyzed by", &self.analyzer)?;
        writeln!(w, " Any additional notes:")?;
        writeln!(w, "    {}", self.notes)?;
        w.flush()?;
        Ok(())
    }

    pub fn read(path: &Path) -> Result<Self> {
        let r = BufReader::new(File::open(path)?);
        let mut d = InfData { band: "Radio".into(), ..Default::default() };
        let mut in_notes = false;
        for l in r.lines() {
            let l = l?;
            if in_notes {
                d.notes.push_str(l.trim());
                continue;
            }
            if l.trim_start().starts_with("Any additional notes") {
                in_notes = true;
                continue;
            }
            let Some((k, v)) = l.split_once('=') else { continue };
            let k = k.trim();
            let v = v.trim();
            match k {
                _ if k.starts_with("Data file name") => d.basename = v.into(),
                _ if k.starts_with("Telescope") => d.telescope = v.into(),
                _ if k.starts_with("Instrument") => d.instrument = v.into(),
                _ if k.starts_with("Object") => d.object = v.into(),
                _ if k.starts_with("J2000 Right") => d.ra = v.into(),
                _ if k.starts_with("J2000 Decl") => d.dec = v.into(),
                _ if k.starts_with("Data observed") => d.observer = v.into(),
                _ if k.starts_with("Epoch") => {
                    let (i, f) = v.split_once('.').unwrap_or((v, "0"));
                    let frac: f64 = format!("0.{f}").parse().unwrap_or(0.0);
                    d.epoch = Epoch::new(i.trim().parse().unwrap_or(0), frac * 86400.0);
                }
                _ if k.starts_with("Barycentered") => d.barycentric = v.starts_with('1'),
                _ if k.starts_with("Number of bins") => d.n = v.parse::<f64>().unwrap_or(0.0) as u64,
                _ if k.starts_with("Width of each") => d.dt = v.parse().unwrap_or(0.0),
                _ if k.starts_with("On/Off") => {
                    let parts: Vec<u64> = v.split(',').filter_map(|x| x.trim().parse().ok()).collect();
                    if parts.len() == 2 {
                        d.breaks.push((parts[0], parts[1]));
                    }
                }
                _ if k.starts_with("Type of observation") => d.band = v.into(),
                _ if k.starts_with("Beam diameter") => d.beam_diam = v.parse().unwrap_or(0.0),
                _ if k.starts_with("Dispersion measure") => d.dm = v.parse().unwrap_or(0.0),
                _ if k.starts_with("Central freq of low") => d.lofreq = v.parse().unwrap_or(0.0),
                _ if k.starts_with("Total bandwidth") => d.bw = v.parse().unwrap_or(0.0),
                _ if k.starts_with("Number of channels") => d.numchan = v.parse().unwrap_or(1),
                _ if k.starts_with("Channel bandwidth") => d.chan_wid = v.parse().unwrap_or(0.0),
                _ if k.starts_with("Data analyzed") => d.analyzer = v.into(),
                _ => {}
            }
        }
        Ok(d)
    }

    pub fn from_header(h: &ObsHeader, basename: &str, dm: f64, n: u64) -> Self {
        let ra = h.ra_hms().map(|(a, b, c)| format!("{a:02}:{b:02}:{c:07.4}")).unwrap_or_else(|| "00:00:00.0000".into());
        let dec = h
            .dec_dms()
            .map(|(neg, a, b, c)| format!("{}{a:02}:{b:02}:{c:07.4}", if neg { "-" } else { "" }))
            .unwrap_or_else(|| "00:00:00.0000".into());
        let lo = h.freqs.mhz.iter().cloned().fold(f64::INFINITY, f64::min);
        InfData {
            basename: basename.into(),
            telescope: h.telescope.clone(),
            instrument: h.backend.clone(),
            object: h.source_name.clone(),
            ra,
            dec,
            observer: "unknown".into(),
            epoch: h.start,
            barycentric: h.barycentric,
            n,
            dt: h.tsamp,
            breaks: vec![],
            band: "Radio".into(),
            beam_diam: 0.0,
            dm,
            lofreq: lo,
            bw: h.freqs.bandwidth(),
            numchan: h.nchan(),
            chan_wid: h.freqs.chan_bw.abs(),
            analyzer: "filterbank".into(),
            notes: h.provenance.join("; "),
        }
    }

    pub fn to_header(&self) -> ObsHeader {
        let freqs = FreqAxis::uniform(self.lofreq, self.chan_wid.max(1e-9), self.numchan.max(1));
        ObsHeader {
            source_name: self.object.clone(),
            telescope: self.telescope.clone(),
            backend: self.instrument.clone(),
            ra_deg: fb_defs::header::ra_from_hms(&self.ra),
            dec_deg: fb_defs::header::dec_from_dms(&self.dec),
            az_deg: None,
            za_deg: None,
            start: self.epoch,
            tsamp: self.dt,
            nsamples: self.n,
            freqs,
            npol: 1,
            pol_type: PolType::Intensity,
            nbits: 32,
            scaling: Scaling::default(),
            barycentric: self.barycentric,
            beam: 0,
            nbeams: 1,
            order: DataOrder::TimePolFreq,
            provenance: vec![format!("presto-inf:{}", self.basename)],
        }
    }
}

pub fn write_dat(root: &Path, s: &Series) -> Result<()> {
    let root = root_of(root);
    let base = root.file_name().and_then(|x| x.to_str()).unwrap_or("data").to_string();
    let mut inf = InfData::from_header(&s.header, &base, s.dm, s.len() as u64);
    inf.lofreq = s.header.freqs.mhz.iter().cloned().fold(f64::INFINITY, f64::min);
    inf.write(&root.with_extension("inf"))?;
    let mut w = BufWriter::new(File::create(root.with_extension("dat"))?);
    for v in &s.data {
        w.write_all(&v.to_le_bytes())?;
    }
    w.flush()?;
    Ok(())
}

pub fn read_dat(root: &Path) -> Result<Series> {
    let root = root_of(root);
    let inf = InfData::read(&root.with_extension("inf"))?;
    let mut bytes = Vec::new();
    File::open(root.with_extension("dat"))?.read_to_end(&mut bytes)?;
    let data: Vec<f32> = bytes.chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect();
    let mut header = inf.to_header();
    header.nsamples = data.len() as u64;
    let f_ref = inf.lofreq + inf.chan_wid * (inf.numchan.max(1) as f64 - 1.0);
    Ok(Series { header, dm: inf.dm, f_ref, data })
}

#[derive(Debug, Clone)]
pub struct PrestoMask {
    pub timesigma: f64,
    pub freqsigma: f64,
    pub mjd: f64,
    pub dtint: f64,
    pub lofreq: f64,
    pub dfreq: f64,
    pub numchan: usize,
    pub numint: usize,
    pub ptsperint: usize,
    pub zap_chans: Vec<i32>,
    pub zap_ints: Vec<i32>,
    pub chans: Vec<Vec<i32>>,
}

fn rd_f64(r: &mut impl Read) -> Result<f64> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    Ok(f64::from_le_bytes(b))
}

fn rd_i32(r: &mut impl Read) -> Result<i32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(i32::from_le_bytes(b))
}

fn rd_vec(r: &mut impl Read, n: usize) -> Result<Vec<i32>> {
    (0..n).map(|_| rd_i32(r)).collect()
}

impl PrestoMask {
    pub fn read(path: &Path) -> Result<Self> {
        let mut r = BufReader::new(File::open(path)?);
        let timesigma = rd_f64(&mut r)?;
        let freqsigma = rd_f64(&mut r)?;
        let mjd = rd_f64(&mut r)?;
        let dtint = rd_f64(&mut r)?;
        let lofreq = rd_f64(&mut r)?;
        let dfreq = rd_f64(&mut r)?;
        let numchan = rd_i32(&mut r)? as usize;
        let numint = rd_i32(&mut r)? as usize;
        let ptsperint = rd_i32(&mut r)? as usize;
        let nzc = rd_i32(&mut r)? as usize;
        let zap_chans = if nzc > 0 { rd_vec(&mut r, nzc)? } else { vec![] };
        let nzi = rd_i32(&mut r)? as usize;
        let zap_ints = if nzi > 0 { rd_vec(&mut r, nzi)? } else { vec![] };
        let counts = rd_vec(&mut r, numint)?;
        let mut chans = Vec::with_capacity(numint);
        for &c in &counts {
            let c = c.max(0) as usize;
            if c > 0 && c < numchan {
                chans.push(rd_vec(&mut r, c)?);
            } else if c >= numchan {
                chans.push((0..numchan as i32).collect());
            } else {
                chans.push(vec![]);
            }
        }
        Ok(Self { timesigma, freqsigma, mjd, dtint, lofreq, dfreq, numchan, numint, ptsperint, zap_chans, zap_ints, chans })
    }

    pub fn write(&self, path: &Path) -> Result<()> {
        let mut w = BufWriter::new(File::create(path)?);
        for v in [self.timesigma, self.freqsigma, self.mjd, self.dtint, self.lofreq, self.dfreq] {
            w.write_all(&v.to_le_bytes())?;
        }
        for v in [self.numchan as i32, self.numint as i32, self.ptsperint as i32] {
            w.write_all(&v.to_le_bytes())?;
        }
        w.write_all(&(self.zap_chans.len() as i32).to_le_bytes())?;
        for v in &self.zap_chans {
            w.write_all(&v.to_le_bytes())?;
        }
        w.write_all(&(self.zap_ints.len() as i32).to_le_bytes())?;
        for v in &self.zap_ints {
            w.write_all(&v.to_le_bytes())?;
        }
        for c in &self.chans {
            w.write_all(&(c.len() as i32).to_le_bytes())?;
        }
        for c in &self.chans {
            if !c.is_empty() && c.len() < self.numchan {
                for v in c {
                    w.write_all(&v.to_le_bytes())?;
                }
            }
        }
        w.flush()?;
        Ok(())
    }

    pub fn to_mask(&self, obs_id: &str) -> RfiMask {
        self.to_mask_with_start(obs_id, None)
    }

    pub fn to_mask_for(&self, h: &ObsHeader) -> RfiMask {
        self.to_mask_with_start(&h.observation_id(), Some((h.start, h.tsamp)))
    }

    fn to_mask_with_start(&self, obs_id: &str, exact: Option<(Epoch, f64)>) -> RfiMask {
        let freqs = FreqAxis::uniform(self.lofreq, self.dfreq, self.numchan);
        let mut start = Epoch::from_mjd(self.mjd);
        if let Some((e, tsamp)) = exact {
            if start.seconds_since(&e).abs() < 0.5 * tsamp {
                start = e;
            }
        }
        let mut m = RfiMask::from_freqs(obs_id.into(), start, self.dtint, self.numint, &freqs, "presto-rfifind");
        for &c in &self.zap_chans {
            if (c as usize) < self.numchan {
                m.zap_channel(c as usize);
            }
        }
        for &i in &self.zap_ints {
            if (i as usize) < self.numint {
                m.zap_interval(i as usize);
            }
        }
        for (i, cs) in self.chans.iter().enumerate() {
            for &c in cs {
                if (c as usize) < self.numchan {
                    m.set(i, c as usize, true);
                }
            }
        }
        m
    }

    pub fn from_mask(m: &RfiMask, tsamp: f64) -> Result<Self> {
        let n = m.nchan();
        let mut order: Vec<usize> = (0..n).collect();
        let centre = |c: usize| 0.5 * (m.chan_lo[c] + m.chan_hi[c]);
        order.sort_by(|a, b| centre(*a).partial_cmp(&centre(*b)).unwrap());
        let dfreq = if n > 1 { centre(order[1]) - centre(order[0]) } else { m.chan_hi[0] - m.chan_lo[0] };
        let lofreq = centre(order[0]);
        for (k, &c) in order.iter().enumerate() {
            if (centre(c) - (lofreq + dfreq * k as f64)).abs() > dfreq.abs() * 1e-3 {
                return Err(FbError::Unsupported("PRESTO masks need uniformly spaced channels".into()));
            }
        }
        let mut chans = Vec::with_capacity(m.nint);
        for i in 0..m.nint {
            chans.push((0..n).filter(|&k| m.get(i, order[k])).map(|k| k as i32).collect::<Vec<i32>>());
        }
        let zap_chans: Vec<i32> = (0..n).filter(|&k| (0..m.nint).all(|i| m.get(i, order[k]))).map(|k| k as i32).collect();
        let zap_ints: Vec<i32> = (0..m.nint).filter(|&i| (0..n).all(|k| m.get(i, order[k]))).map(|i| i as i32).collect();
        Ok(Self {
            timesigma: 10.0,
            freqsigma: 4.0,
            mjd: m.start.mjd(),
            dtint: m.interval,
            lofreq,
            dfreq,
            numchan: n,
            numint: m.nint,
            ptsperint: (m.interval / tsamp).round().max(1.0) as usize,
            zap_chans,
            zap_ints,
            chans,
        })
    }
}
