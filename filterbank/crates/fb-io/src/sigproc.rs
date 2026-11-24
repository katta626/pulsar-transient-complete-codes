use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::Arc;

use fb_defs::{Block, DataOrder, Epoch, FbError, FreqAxis, ObsHeader, PolType, Result, Scaling, Series};
use memmap2::Mmap;

use crate::bits::{pack, unpack_into, BitOrder};
use crate::reader::Reader;

const TELESCOPES: &[(i32, &str)] = &[
    (0, "Fake"),
    (1, "Arecibo"),
    (2, "Ooty"),
    (3, "Nancay"),
    (4, "Parkes"),
    (5, "Jodrell"),
    (6, "GBT"),
    (7, "GMRT"),
    (8, "Effelsberg"),
    (9, "140ft"),
    (10, "SRT"),
    (64, "MeerKAT"),
    (65, "KAT-7"),
    (82, "eMerlin"),
];

const MACHINES: &[(i32, &str)] = &[
    (0, "FAKE"),
    (1, "PSPM"),
    (2, "WAPP"),
    (3, "AOFTM"),
    (4, "BPP"),
    (5, "OOTY"),
    (6, "SCAMP"),
    (7, "GMRTFB"),
    (8, "PULSAR2000"),
    (9, "PARSPEC"),
    (10, "BPSR"),
    (14, "GMRTNEW"),
    (64, "KAT"),
    (65, "KAT-DC2"),
    (82, "loft-e"),
    (83, "ROACH"),
    (88, "AT-PAF"),
];

pub fn telescope_name(id: i32) -> String {
    TELESCOPES.iter().find(|t| t.0 == id).map(|t| t.1.to_string()).unwrap_or_else(|| format!("telescope_{id}"))
}

pub fn telescope_id(name: &str) -> i32 {
    let n = name.trim().to_ascii_lowercase();
    if let Some(rest) = n.strip_prefix("telescope_") {
        return rest.parse().unwrap_or(0);
    }
    TELESCOPES.iter().find(|t| t.1.to_ascii_lowercase() == n || n.contains(&t.1.to_ascii_lowercase())).map(|t| t.0).unwrap_or(0)
}

pub fn machine_name(id: i32) -> String {
    MACHINES.iter().find(|t| t.0 == id).map(|t| t.1.to_string()).unwrap_or_else(|| format!("machine_{id}"))
}

pub fn machine_id(name: &str) -> i32 {
    let n = name.trim().to_ascii_lowercase();
    if let Some(rest) = n.strip_prefix("machine_") {
        return rest.parse().unwrap_or(0);
    }
    MACHINES.iter().find(|t| t.1.to_ascii_lowercase() == n).map(|t| t.0).unwrap_or(0)
}

fn packed_to_deg(v: f64, hours: bool) -> f64 {
    let neg = v < 0.0;
    let a = v.abs();
    let d = (a / 10000.0).floor();
    let m = ((a - d * 10000.0) / 100.0).floor();
    let s = a - d * 10000.0 - m * 100.0;
    let x = d + m / 60.0 + s / 3600.0;
    let x = if hours { x * 15.0 } else { x };
    if neg { -x } else { x }
}

fn deg_to_packed(v: f64, hours: bool) -> f64 {
    let neg = v < 0.0;
    let a = if hours { v.rem_euclid(360.0) / 15.0 } else { v.abs() };
    let d = a.floor();
    let mf = (a - d) * 60.0;
    let m = mf.floor();
    let s = (mf - m) * 60.0;
    let p = d * 10000.0 + m * 100.0 + s;
    if neg { -p } else { p }
}

#[derive(Debug, Clone, Default)]
pub struct SigprocHeader {
    pub ints: Vec<(String, i64)>,
    pub doubles: Vec<(String, f64)>,
    pub strings: Vec<(String, String)>,
    pub freq_table: Vec<f64>,
    pub signed: Option<bool>,
    pub header_len: usize,
}

impl SigprocHeader {
    pub fn int(&self, k: &str) -> Option<i64> {
        self.ints.iter().find(|x| x.0 == k).map(|x| x.1)
    }
    pub fn double(&self, k: &str) -> Option<f64> {
        self.doubles.iter().find(|x| x.0 == k).map(|x| x.1)
    }
    pub fn string(&self, k: &str) -> Option<String> {
        self.strings.iter().find(|x| x.0 == k).map(|x| x.1.clone())
    }
}

fn read_str(buf: &[u8], pos: &mut usize) -> Result<String> {
    if *pos + 4 > buf.len() {
        return Err(FbError::Format("truncated SIGPROC header".into()));
    }
    let n = i32::from_le_bytes(buf[*pos..*pos + 4].try_into().unwrap());
    if !(1..=80).contains(&n) {
        return Err(FbError::Format(format!("bad SIGPROC string length {n}")));
    }
    *pos += 4;
    let s = String::from_utf8_lossy(&buf[*pos..*pos + n as usize]).to_string();
    *pos += n as usize;
    Ok(s)
}

pub fn parse_header(buf: &[u8]) -> Result<SigprocHeader> {
    let mut pos = 0usize;
    if read_str(buf, &mut pos)? != "HEADER_START" {
        return Err(FbError::Format("missing HEADER_START".into()));
    }
    let mut h = SigprocHeader::default();
    let take = |pos: &mut usize, n: usize| -> Result<&[u8]> {
        if *pos + n > buf.len() {
            return Err(FbError::Format("truncated SIGPROC header".into()));
        }
        let s = &buf[*pos..*pos + n];
        *pos += n;
        Ok(s)
    };
    loop {
        let key = read_str(buf, &mut pos)?;
        match key.as_str() {
            "HEADER_END" => break,
            "rawdatafile" | "source_name" => {
                let v = read_str(buf, &mut pos)?;
                h.strings.push((key, v));
            }
            "FREQUENCY_START" | "FREQUENCY_END" => {}
            "fchannel" => {
                let b = take(&mut pos, 8)?;
                h.freq_table.push(f64::from_le_bytes(b.try_into().unwrap()));
            }
            "telescope_id" | "machine_id" | "data_type" | "barycentric" | "pulsarcentric" | "nbits" | "nsamples" | "nchans"
            | "nifs" | "nbins" | "nbeams" | "ibeam" => {
                let b = take(&mut pos, 4)?;
                h.ints.push((key, i32::from_le_bytes(b.try_into().unwrap()) as i64));
            }
            "npuls" => {
                let b = take(&mut pos, 8)?;
                h.ints.push((key, i64::from_le_bytes(b.try_into().unwrap())));
            }
            "signed" => {
                let b = take(&mut pos, 1)?;
                h.signed = Some(b[0] != 0);
            }
            "az_start" | "za_start" | "src_raj" | "src_dej" | "tstart" | "tsamp" | "fch1" | "foff" | "refdm" | "period"
            | "refrm" => {
                let b = take(&mut pos, 8)?;
                h.doubles.push((key, f64::from_le_bytes(b.try_into().unwrap())));
            }
            other => return Err(FbError::Format(format!("unknown SIGPROC keyword '{other}'"))),
        }
    }
    h.header_len = pos;
    Ok(h)
}

pub struct SigprocReader {
    mmap: Mmap,
    pub raw: SigprocHeader,
    header: Arc<ObsHeader>,
    nbits: u32,
    signed: bool,
}

impl SigprocReader {
    pub fn open(path: &Path) -> Result<Self> {
        let f = File::open(path)?;
        let mmap = unsafe { Mmap::map(&f)? };
        let raw = parse_header(&mmap)?;
        let nbits = raw.int("nbits").unwrap_or(8) as u32;
        let nchan = raw.int("nchans").unwrap_or(1) as usize;
        let nifs = raw.int("nifs").unwrap_or(1).max(1) as usize;
        let signed = raw.signed.unwrap_or(false);
        let bytes = mmap.len() - raw.header_len;
        let nsamples = (bytes as u64 * 8) / (nbits as u64 * nchan as u64 * nifs as u64);
        let freqs = if raw.freq_table.len() == nchan {
            FreqAxis::from_centres(raw.freq_table.clone(), None)?
        } else {
            FreqAxis::uniform(raw.double("fch1").unwrap_or(0.0), raw.double("foff").unwrap_or(-1.0), nchan)
        };
        let raj = raw.double("src_raj");
        let dej = raw.double("src_dej");
        let header = ObsHeader {
            source_name: raw.string("source_name").unwrap_or_else(|| "UNKNOWN".into()),
            telescope: telescope_name(raw.int("telescope_id").unwrap_or(0) as i32),
            backend: machine_name(raw.int("machine_id").unwrap_or(0) as i32),
            ra_deg: raj.filter(|v| *v != 0.0 && *v != -1.0).map(|v| packed_to_deg(v, true)),
            dec_deg: dej.filter(|v| *v != -1.0).map(|v| packed_to_deg(v, false)),
            az_deg: raw.double("az_start").filter(|v| *v != -1.0),
            za_deg: raw.double("za_start").filter(|v| *v != -1.0),
            start: Epoch::from_mjd(raw.double("tstart").unwrap_or(0.0)),
            tsamp: raw.double("tsamp").unwrap_or(0.0),
            nsamples,
            freqs,
            npol: nifs,
            pol_type: if nifs == 1 { PolType::Intensity } else { PolType::AABB },
            nbits,
            scaling: Scaling { zero_off: 0.0, signed, applied: false },
            barycentric: raw.int("barycentric").unwrap_or(0) == 1,
            beam: raw.int("ibeam").unwrap_or(1).max(0) as u32,
            nbeams: raw.int("nbeams").unwrap_or(1).max(1) as u32,
            order: DataOrder::TimePolFreq,
            provenance: vec![format!("sigproc:{}", path.file_name().and_then(|s| s.to_str()).unwrap_or(""))],
        };
        Ok(Self { mmap, raw, header: Arc::new(header), nbits, signed })
    }
}

impl Reader for SigprocReader {
    fn header(&self) -> Arc<ObsHeader> {
        self.header.clone()
    }

    fn read(&mut self, start: u64, nsamp: usize) -> Result<Block> {
        let h = self.header.clone();
        let end = (start + nsamp as u64).min(h.nsamples);
        let n = end.saturating_sub(start) as usize;
        let per = h.npol * h.nchan();
        let bits_per_samp = per as u64 * self.nbits as u64;
        let b0 = self.raw.header_len + (start * bits_per_samp / 8) as usize;
        let b1 = self.raw.header_len + (end * bits_per_samp / 8) as usize;
        let mut data = Vec::with_capacity(n * per);
        unpack_into(&self.mmap[b0..b1], self.nbits, BitOrder::LsbFirst, self.signed && self.nbits <= 8, &mut data);
        data.truncate(n * per);
        Block::new(h, start, n, data)
    }
}

pub struct SigprocWriter {
    out: BufWriter<File>,
    nbits: u32,
    signed: bool,
    width: usize,
    pub samples_written: u64,
}

fn put_str(buf: &mut Vec<u8>, s: &str) {
    buf.extend((s.len() as i32).to_le_bytes());
    buf.extend(s.as_bytes());
}

fn put_i32(buf: &mut Vec<u8>, k: &str, v: i32) {
    put_str(buf, k);
    buf.extend(v.to_le_bytes());
}

fn put_f64(buf: &mut Vec<u8>, k: &str, v: f64) {
    put_str(buf, k);
    buf.extend(v.to_le_bytes());
}

pub struct WriteOptions {
    pub nbits: u32,
    pub data_type: i32,
    pub refdm: Option<f64>,
    pub rawdatafile: String,
}

pub fn header_bytes(h: &ObsHeader, opts: &WriteOptions) -> Vec<u8> {
    let mut b = Vec::new();
    put_str(&mut b, "HEADER_START");
    put_str(&mut b, "rawdatafile");
    put_str(&mut b, &opts.rawdatafile.chars().take(80).collect::<String>());
    put_str(&mut b, "source_name");
    put_str(&mut b, &h.source_name.chars().take(80).collect::<String>());
    put_i32(&mut b, "machine_id", machine_id(&h.backend));
    put_i32(&mut b, "telescope_id", telescope_id(&h.telescope));
    if let Some(ra) = h.ra_deg {
        put_f64(&mut b, "src_raj", deg_to_packed(ra, true));
    }
    if let Some(dec) = h.dec_deg {
        put_f64(&mut b, "src_dej", deg_to_packed(dec, false));
    }
    if let Some(az) = h.az_deg {
        put_f64(&mut b, "az_start", az);
    }
    if let Some(za) = h.za_deg {
        put_f64(&mut b, "za_start", za);
    }
    put_i32(&mut b, "data_type", opts.data_type);
    if opts.data_type == 2 {
        put_f64(&mut b, "refdm", opts.refdm.unwrap_or(0.0));
        put_f64(&mut b, "fch1", h.freqs.fch1());
        put_i32(&mut b, "nchans", 1);
    } else {
        put_f64(&mut b, "fch1", h.freqs.fch1());
        put_f64(&mut b, "foff", h.freqs.foff());
        put_i32(&mut b, "nchans", h.nchan() as i32);
    }
    put_i32(&mut b, "nbeams", h.nbeams as i32);
    put_i32(&mut b, "ibeam", h.beam as i32);
    put_i32(&mut b, "nbits", opts.nbits as i32);
    put_f64(&mut b, "tstart", h.start.mjd());
    put_f64(&mut b, "tsamp", h.tsamp);
    put_i32(&mut b, "nifs", if opts.data_type == 2 { 1 } else { h.npol as i32 });
    put_i32(&mut b, "barycentric", h.barycentric as i32);
    if opts.nbits == 8 {
        put_str(&mut b, "signed");
        b.push(h.scaling.signed as u8);
    }
    put_str(&mut b, "HEADER_END");
    b
}

impl SigprocWriter {
    pub fn create(path: &Path, h: &ObsHeader, opts: WriteOptions) -> Result<Self> {
        if ![1, 2, 4, 8, 16, 32].contains(&opts.nbits) {
            return Err(FbError::Unsupported(format!("nbits={}", opts.nbits)));
        }
        if !h.freqs.is_uniform(1e-3) {
            return Err(FbError::Unsupported("SIGPROC needs uniformly spaced channels".into()));
        }
        let mut out = BufWriter::with_capacity(1 << 22, File::create(path)?);
        out.write_all(&header_bytes(h, &opts))?;
        let width = if opts.data_type == 2 { 1 } else { h.npol * h.nchan() };
        Ok(Self { out, nbits: opts.nbits, signed: h.scaling.signed, width, samples_written: 0 })
    }

    pub fn write_values(&mut self, values: &[f32]) -> Result<()> {
        if values.len() % self.width != 0 {
            return Err(FbError::Invalid("value count is not a whole number of spectra".into()));
        }
        self.out.write_all(&pack(values, self.nbits, BitOrder::LsbFirst, self.signed))?;
        self.samples_written += (values.len() / self.width) as u64;
        Ok(())
    }

    pub fn write_block(&mut self, b: &Block) -> Result<()> {
        self.write_values(&b.data)
    }

    pub fn finish(mut self) -> Result<u64> {
        self.out.flush()?;
        Ok(self.samples_written)
    }
}

pub fn write_tim(path: &Path, s: &Series) -> Result<()> {
    let opts = WriteOptions { nbits: 32, data_type: 2, refdm: Some(s.dm), rawdatafile: s.header.source_name.clone() };
    let mut h = s.header.clone();
    h.freqs = FreqAxis::uniform(s.f_ref, h.freqs.chan_bw, 1);
    h.npol = 1;
    let mut w = SigprocWriter::create(path, &h, opts)?;
    w.write_values(&s.data)?;
    w.finish()?;
    Ok(())
}
