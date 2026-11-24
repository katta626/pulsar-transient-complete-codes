use std::collections::HashMap;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use fb_defs::{Block, DataOrder, Epoch, FbError, FreqAxis, ObsHeader, PolType, Result, Scaling};
use memmap2::Mmap;

use crate::bits::{unpack_into, BitOrder};
use crate::reader::Reader;

pub fn parse_header(text: &str) -> HashMap<String, String> {
    let mut m = HashMap::new();
    for line in text.lines() {
        let l = line.split('#').next().unwrap_or("").trim();
        if l.is_empty() {
            continue;
        }
        let mut it = l.splitn(2, char::is_whitespace);
        if let (Some(k), Some(v)) = (it.next(), it.next()) {
            m.insert(k.trim().to_ascii_uppercase(), v.trim().to_string());
        }
    }
    m
}

pub fn civil_to_mjd(y: i64, mo: i64, d: i64) -> i64 {
    let a = (14 - mo) / 12;
    let yy = y + 4800 - a;
    let mm = mo + 12 * a - 3;
    let jdn = d + (153 * mm + 2) / 5 + 365 * yy + yy / 4 - yy / 100 + yy / 400 - 32045;
    jdn - 2_400_001
}

pub fn utc_to_epoch(s: &str) -> Option<Epoch> {
    let s = s.trim().replace('T', "-").replace(' ', "-");
    let parts: Vec<&str> = s.split(['-', ':']).collect();
    if parts.len() < 6 {
        return None;
    }
    let y = parts[0].parse().ok()?;
    let mo = parts[1].parse().ok()?;
    let d = parts[2].parse().ok()?;
    let hh: f64 = parts[3].parse().ok()?;
    let mi: f64 = parts[4].parse().ok()?;
    let ss: f64 = parts[5].parse().ok()?;
    Some(Epoch::new(civil_to_mjd(y, mo, d), hh * 3600.0 + mi * 60.0 + ss))
}

pub struct DadaReader {
    mmap: Mmap,
    hdr_size: usize,
    pub keys: HashMap<String, String>,
    header: Arc<ObsHeader>,
    nbits: u32,
    order: BitOrder,
    npol_raw: usize,
    sum_pols: bool,
}

impl DadaReader {
    pub fn open(path: &Path, sum_pols: bool) -> Result<Self> {
        let f = File::open(path)?;
        let mmap = unsafe { Mmap::map(&f)? };
        let probe = String::from_utf8_lossy(&mmap[..mmap.len().min(65536)]).to_string();
        let first = parse_header(&probe);
        let hdr_size: usize = first.get("HDR_SIZE").and_then(|v| v.parse().ok()).unwrap_or(4096);
        let text = String::from_utf8_lossy(&mmap[..hdr_size.min(mmap.len())]).trim_end_matches('\0').to_string();
        let k = parse_header(&text);
        let get = |name: &str| k.get(name).cloned();
        let num = |name: &str| get(name).and_then(|v| v.parse::<f64>().ok());
        let ndim = num("NDIM").unwrap_or(1.0) as usize;
        if ndim != 1 {
            return Err(FbError::Unsupported("DADA with NDIM=2 is voltage data; channelise and detect it first (e.g. with dspsr)".into()));
        }
        let nchan = num("NCHAN").ok_or_else(|| FbError::Format("DADA header has no NCHAN".into()))? as usize;
        let nbits = num("NBIT").unwrap_or(8.0) as u32;
        let npol_raw = num("NPOL").unwrap_or(1.0) as usize;
        let tsamp = num("TSAMP").ok_or_else(|| FbError::Format("DADA header has no TSAMP".into()))? * 1e-6;
        let bw = num("BW").or_else(|| num("BANDWIDTH")).unwrap_or(-(nchan as f64));
        let fc = num("FREQ").or_else(|| num("CFREQ")).unwrap_or(1400.0);
        let chan_bw = bw / nchan as f64;
        let fch1 = fc - bw / 2.0 + chan_bw / 2.0;
        let bytes_per_samp = nchan * npol_raw * nbits as usize / 8;
        let obs_offset = num("OBS_OFFSET").unwrap_or(0.0) as u64;
        let utc = get("UTC_START").and_then(|s| utc_to_epoch(&s)).ok_or_else(|| FbError::Format("DADA header has no usable UTC_START".into()))?;
        let start = utc.add_seconds((obs_offset / bytes_per_samp.max(1) as u64) as f64 * tsamp);
        let data_bytes = mmap.len().saturating_sub(hdr_size);
        let nsamples = (data_bytes / bytes_per_samp.max(1)) as u64;
        let order = match get("BIT_ORDER").map(|s| s.to_ascii_uppercase()) {
            Some(s) if s.starts_with("MSB") => BitOrder::MsbFirst,
            _ => BitOrder::LsbFirst,
        };
        let header = ObsHeader {
            source_name: get("SOURCE").unwrap_or_else(|| "UNKNOWN".into()),
            telescope: get("TELESCOPE").unwrap_or_default(),
            backend: get("INSTRUMENT").or_else(|| get("MACHINE")).unwrap_or_else(|| "DADA".into()),
            ra_deg: get("RA").and_then(|s| fb_defs::header::ra_from_hms(&s)),
            dec_deg: get("DEC").and_then(|s| fb_defs::header::dec_from_dms(&s)),
            az_deg: None,
            za_deg: None,
            start,
            tsamp,
            nsamples,
            freqs: FreqAxis::uniform(fch1, chan_bw, nchan),
            npol: if sum_pols { 1 } else { npol_raw },
            pol_type: if npol_raw == 1 || sum_pols { PolType::Intensity } else { PolType::AABB },
            nbits,
            scaling: Scaling::default(),
            barycentric: false,
            beam: num("BEAM").or_else(|| num("IBEAM")).unwrap_or(0.0) as u32,
            nbeams: num("NBEAM").unwrap_or(1.0) as u32,
            order: DataOrder::TimePolFreq,
            provenance: vec![format!("dada:{}", path.file_name().and_then(|s| s.to_str()).unwrap_or(""))],
        };
        Ok(Self { mmap, hdr_size, keys: k, header: Arc::new(header), nbits, order, npol_raw, sum_pols })
    }
}

impl Reader for DadaReader {
    fn header(&self) -> Arc<ObsHeader> {
        self.header.clone()
    }

    fn read(&mut self, start: u64, nsamp: usize) -> Result<Block> {
        let h = self.header.clone();
        let end = (start + nsamp as u64).min(h.nsamples);
        let n = end.saturating_sub(start) as usize;
        let nchan = h.nchan();
        let raw_w = nchan * self.npol_raw;
        let bps = raw_w * self.nbits as usize / 8;
        let b0 = self.hdr_size + start as usize * bps;
        let mut vals = Vec::with_capacity(n * raw_w);
        unpack_into(&self.mmap[b0..b0 + n * bps], self.nbits, self.order, false, &mut vals);
        vals.truncate(n * raw_w);
        if self.sum_pols && self.npol_raw > 1 {
            let mut out = vec![0f32; n * nchan];
            for t in 0..n {
                for p in 0..self.npol_raw.min(2) {
                    for c in 0..nchan {
                        out[t * nchan + c] += vals[(t * self.npol_raw + p) * nchan + c];
                    }
                }
            }
            return Block::new(h, start, n, out);
        }
        Block::new(h, start, n, vals)
    }
}
