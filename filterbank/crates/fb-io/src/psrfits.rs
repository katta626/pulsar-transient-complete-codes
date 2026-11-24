use std::path::{Path, PathBuf};
use std::sync::Arc;

use fb_defs::header::{dec_from_dms, ra_from_hms};
use fb_defs::{Block, DataOrder, Epoch, FbError, FreqAxis, ObsHeader, PolType, Result, Scaling};

use crate::bits::{unpack_into, BitOrder};
use crate::fits::{FitsFile, Table};
use crate::reader::{Gap, Reader};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ZeroOffPolicy {
    Spec,
    Dspsr,
    Ignore,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScalePolicy {
    Auto,
    Always,
    Never,
}

#[derive(Debug, Clone, Copy)]
pub struct PsrfitsOptions {
    pub zero_off: ZeroOffPolicy,
    pub scales: ScalePolicy,
    pub weights: bool,
    pub sum_pols: bool,
}

impl Default for PsrfitsOptions {
    fn default() -> Self {
        Self { zero_off: ZeroOffPolicy::Spec, scales: ScalePolicy::Auto, weights: true, sum_pols: true }
    }
}

impl PsrfitsOptions {
    pub fn raw() -> Self {
        Self { zero_off: ZeroOffPolicy::Ignore, scales: ScalePolicy::Never, weights: false, sum_pols: false }
    }
}

struct RowRef {
    file: usize,
    row: usize,
    start: u64,
}

pub struct PsrfitsReader {
    paths: Vec<PathBuf>,
    files: Vec<FitsFile>,
    rows: Vec<RowRef>,
    header: Arc<ObsHeader>,
    opts: PsrfitsOptions,
    nsblk: usize,
    nbits: u32,
    npol_raw: usize,
    nchan: usize,
    signed: bool,
    zero_off: f64,
    pol_type: PolType,
    use_scales: bool,
    gaps: Vec<Gap>,
    parkes_1bit_flip: bool,
}

fn primary_epoch(f: &FitsFile) -> Result<Epoch> {
    let p = f.primary();
    Ok(Epoch::from_parts(p.req_i64("STT_IMJD")?, p.f64("STT_SMJD").unwrap_or(0.0), p.f64("STT_OFFS").unwrap_or(0.0)))
}

impl PsrfitsReader {
    pub fn open<P: AsRef<Path>>(paths: &[P], opts: PsrfitsOptions) -> Result<Self> {
        if paths.is_empty() {
            return Err(FbError::Invalid("no PSRFITS files given".into()));
        }
        let mut items: Vec<(PathBuf, FitsFile, Epoch)> = Vec::new();
        for p in paths {
            let f = FitsFile::open(p.as_ref())?;
            let mode = f.primary().str("OBS_MODE").unwrap_or_default().to_ascii_uppercase();
            if !mode.starts_with("SEARCH") {
                return Err(FbError::Unsupported(format!(
                    "{} has OBS_MODE '{}'; only SEARCH mode is read as filterbank data",
                    p.as_ref().display(),
                    mode
                )));
            }
            let e = primary_epoch(&f)?;
            items.push((p.as_ref().to_path_buf(), f, e));
        }
        items.sort_by(|a, b| a.2.seconds_since(&b.2).partial_cmp(&0.0).unwrap());
        let (paths, files, epochs): (Vec<_>, Vec<_>, Vec<_>) = items.into_iter().fold((vec![], vec![], vec![]), |mut acc, (p, f, e)| {
            acc.0.push(p);
            acc.1.push(f);
            acc.2.push(e);
            acc
        });

        let first = &files[0];
        let prim = first.primary();
        let sub = first.table("SUBINT")?;
        let sh = sub.header;
        let nbits = sh.req_i64("NBITS")? as u32;
        let nchan = sh.req_i64("NCHAN")? as usize;
        let npol_raw = sh.req_i64("NPOL")? as usize;
        let nsblk = sh.req_i64("NSBLK")? as usize;
        let tbin = sh.req_f64("TBIN")?;
        let signed = sh.i64("SIGNINT").unwrap_or(0) == 1;
        let zero_off = sh.f64("ZERO_OFF").unwrap_or(0.0).abs();
        let pol_type = PolType::parse(&sh.str("POL_TYPE").unwrap_or_default(), npol_raw);
        if ![1, 2, 4, 8, 16, 32].contains(&nbits) {
            return Err(FbError::Unsupported(format!("NBITS={nbits}")));
        }

        let freqs = match sub.f64s(0, "DAT_FREQ") {
            Some(f) if f.len() == nchan => FreqAxis::from_centres(f, sh.f64("CHAN_BW"))?,
            _ => {
                let bw = sh.f64("CHAN_BW").or_else(|| prim.f64("OBSBW").map(|b| b / nchan as f64)).unwrap_or(-1.0);
                let fc = prim.f64("OBSFREQ").unwrap_or(0.0);
                let half = (nchan / 2) as f64;
                FreqAxis::uniform(fc - bw * half, bw, nchan)
            }
        };

        let mut rows = Vec::new();
        let start0 = epochs[0];
        for (fi, f) in files.iter().enumerate() {
            let t = f.table("SUBINT")?;
            let h = t.header;
            for (k, expect) in [("NBITS", nbits as i64), ("NCHAN", nchan as i64), ("NPOL", npol_raw as i64), ("NSBLK", nsblk as i64)] {
                let got = h.req_i64(k)?;
                if got != expect {
                    return Err(FbError::Incompatible(format!("{} has {k}={got}, expected {expect}", paths[fi].display())));
                }
            }
            if let Some(fr) = t.f64s(0, "DAT_FREQ") {
                let other = FreqAxis::from_centres(fr, Some(freqs.chan_bw))?;
                if !freqs.same_as(&other, 1e-3) {
                    return Err(FbError::Incompatible(format!("{} has a different DAT_FREQ axis", paths[fi].display())));
                }
            }
            let file_off = epochs[fi].seconds_since(&start0) / tbin;
            let nsuboffs = h.i64("NSUBOFFS").unwrap_or(0).max(0) as f64;
            let offs: Vec<f64> = (0..t.nrows).map(|r| t.f64_scalar(r, "OFFS_SUB").unwrap_or(0.0)).collect();
            let all_zero = offs.iter().all(|o| *o == 0.0);
            let tsub = (nsblk as f64) * tbin;
            for (r, o) in offs.iter().enumerate() {
                let local = if all_zero { (nsuboffs + r as f64) * nsblk as f64 } else { (o - tsub / 2.0) / tbin };
                let s = (file_off + local).round().max(0.0) as u64;
                rows.push(RowRef { file: fi, row: r, start: s });
            }
        }
        rows.sort_by_key(|r| r.start);
        let mut gaps = Vec::new();
        for w in rows.windows(2) {
            let end = w[0].start + nsblk as u64;
            if w[1].start > end {
                gaps.push(Gap { start: end, len: w[1].start - end });
            } else if w[1].start < end {
                return Err(FbError::Format(format!("overlapping SUBINT rows near sample {}", w[1].start)));
            }
        }
        let nsamples = rows.last().map(|r| r.start + nsblk as u64).unwrap_or(0);

        let use_scales = match opts.scales {
            ScalePolicy::Never => false,
            ScalePolicy::Always => true,
            ScalePolicy::Auto => scales_nontrivial(&sub, opts.weights),
        };

        let telescope = prim.str("TELESCOP").unwrap_or_default();
        let backend = prim.str("BACKEND").unwrap_or_default();
        let parkes_1bit_flip = nbits == 1 && telescope.to_ascii_uppercase().contains("PARKES") && backend.to_ascii_uppercase() == "FB_1BIT";
        let npol_out = if opts.sum_pols { 1 } else { npol_raw };
        let applied = use_scales || opts.zero_off != ZeroOffPolicy::Ignore;
        let header = ObsHeader {
            source_name: prim.str("SRC_NAME").unwrap_or_else(|| "UNKNOWN".into()),
            telescope,
            backend,
            ra_deg: prim.str("RA").and_then(|s| ra_from_hms(&s)),
            dec_deg: prim.str("DEC").and_then(|s| dec_from_dms(&s)),
            az_deg: sub.f64_scalar(0, "TEL_AZ"),
            za_deg: sub.f64_scalar(0, "TEL_ZEN"),
            start: start0,
            tsamp: tbin,
            nsamples,
            freqs,
            npol: npol_out,
            pol_type: if opts.sum_pols { PolType::Intensity } else { pol_type },
            nbits,
            scaling: Scaling { zero_off, signed, applied },
            barycentric: false,
            beam: prim.i64("IBEAM").unwrap_or(1).max(0) as u32,
            nbeams: prim.i64("NBEAMS").unwrap_or(1).max(1) as u32,
            order: DataOrder::TimePolFreq,
            provenance: vec![format!(
                "psrfits:{} files={} zero_off={:?} scales={} weights={} sum_pols={} gaps={}",
                paths[0].file_name().and_then(|s| s.to_str()).unwrap_or(""),
                paths.len(),
                opts.zero_off,
                use_scales,
                opts.weights && use_scales,
                opts.sum_pols,
                gaps.len()
            )],
        };
        Ok(Self {
            paths,
            files,
            rows,
            header: Arc::new(header),
            opts,
            nsblk,
            nbits,
            npol_raw,
            nchan,
            signed,
            zero_off,
            pol_type,
            use_scales,
            gaps,
            parkes_1bit_flip,
        })
    }

    pub fn paths(&self) -> &[PathBuf] {
        &self.paths
    }

    pub fn nsblk(&self) -> usize {
        self.nsblk
    }

    fn zero_for_bits(&self) -> f64 {
        match self.opts.zero_off {
            ZeroOffPolicy::Ignore => 0.0,
            ZeroOffPolicy::Spec => self.zero_off,
            ZeroOffPolicy::Dspsr => {
                if self.nbits == 2 || self.nbits == 16 {
                    0.0
                } else {
                    self.zero_off
                }
            }
        }
    }

    fn decode_row(&self, rr: &RowRef) -> Result<Vec<f32>> {
        let t = self.files[rr.file].table("SUBINT")?;
        let col = t.column("DATA").ok_or_else(|| FbError::Format("SUBINT has no DATA column".into()))?;
        let raw = t.cell(rr.row, col);
        let n = self.nsblk * self.npol_raw * self.nchan;
        let mut vals: Vec<f32> = Vec::with_capacity(n);
        match (self.nbits, col.code) {
            (16, 'I') => {
                let tzero = t.header.f64(&format!("TZERO{}", t.columns.iter().position(|c| c.name == "DATA").unwrap() + 1)).unwrap_or(0.0);
                vals.extend(raw.chunks_exact(2).map(|c| i16::from_be_bytes([c[0], c[1]]) as f32 + tzero as f32));
            }
            (16, _) => vals.extend(raw.chunks_exact(2).map(|c| {
                if self.signed { i16::from_be_bytes([c[0], c[1]]) as f32 } else { u16::from_be_bytes([c[0], c[1]]) as f32 }
            })),
            (32, _) => vals.extend(raw.chunks_exact(4).map(|c| f32::from_be_bytes([c[0], c[1], c[2], c[3]]))),
            (1, _) if self.parkes_1bit_flip => {
                let mut flipped = raw.to_vec();
                for chunk in flipped.chunks_mut(8) {
                    chunk.reverse();
                }
                unpack_into(&flipped, 1, BitOrder::MsbFirst, false, &mut vals);
            }
            (nb, _) => unpack_into(raw, nb, BitOrder::MsbFirst, self.signed && nb == 8, &mut vals),
        }
        vals.truncate(n);
        if vals.len() < n {
            return Err(FbError::Format(format!("row {} of file {} is short", rr.row, rr.file)));
        }

        let zero = self.zero_for_bits() as f32;
        if self.use_scales {
            let scl = t.f64s(rr.row, "DAT_SCL").unwrap_or_else(|| vec![1.0; self.nchan * self.npol_raw]);
            let offs = t.f64s(rr.row, "DAT_OFFS").unwrap_or_else(|| vec![0.0; self.nchan * self.npol_raw]);
            let wts = if self.opts.weights { t.f64s(rr.row, "DAT_WTS").unwrap_or_else(|| vec![1.0; self.nchan]) } else { vec![1.0; self.nchan] };
            for s in 0..self.nsblk {
                for p in 0..self.npol_raw {
                    for c in 0..self.nchan {
                        let i = (s * self.npol_raw + p) * self.nchan + c;
                        let k = (p * self.nchan + c).min(scl.len().saturating_sub(1));
                        let ko = (p * self.nchan + c).min(offs.len().saturating_sub(1));
                        vals[i] = (((vals[i] - zero) as f64 * scl[k] + offs[ko]) * wts[c.min(wts.len() - 1)]) as f32;
                    }
                }
            }
        } else if zero != 0.0 {
            for v in vals.iter_mut() {
                *v -= zero;
            }
        }

        if self.opts.sum_pols && self.npol_raw > 1 {
            let terms = self.pol_type.total_intensity_terms(self.npol_raw);
            let mut out = vec![0f32; self.nsblk * self.nchan];
            for s in 0..self.nsblk {
                for &p in &terms {
                    let base = (s * self.npol_raw + p) * self.nchan;
                    for c in 0..self.nchan {
                        out[s * self.nchan + c] += vals[base + c];
                    }
                }
            }
            return Ok(out);
        }
        Ok(vals)
    }

    fn gap_fill(&self) -> f32 {
        if self.opts.zero_off == ZeroOffPolicy::Ignore && !self.use_scales {
            (self.zero_off.round() as f32) * if self.opts.sum_pols && self.npol_raw > 1 { 2.0 } else { 1.0 }
        } else {
            0.0
        }
    }
}

fn scales_nontrivial(t: &Table<'_>, weights: bool) -> bool {
    let check = |name: &str, unit: f64| t.f64s(0, name).map(|v| v.iter().any(|x| (x - unit).abs() > 1e-9)).unwrap_or(false);
    check("DAT_SCL", 1.0) || check("DAT_OFFS", 0.0) || (weights && check("DAT_WTS", 1.0))
}

impl Reader for PsrfitsReader {
    fn header(&self) -> Arc<ObsHeader> {
        self.header.clone()
    }

    fn gaps(&self) -> &[Gap] {
        &self.gaps
    }

    fn read(&mut self, start: u64, nsamp: usize) -> Result<Block> {
        let h = self.header.clone();
        let end = (start + nsamp as u64).min(h.nsamples);
        let n = end.saturating_sub(start) as usize;
        let width = h.npol * h.nchan();
        let mut data = vec![self.gap_fill(); n * width];
        let first = self.rows.partition_point(|r| r.start + self.nsblk as u64 <= start);
        for rr in &self.rows[first..] {
            if rr.start >= end {
                break;
            }
            let vals = self.decode_row(rr)?;
            let r0 = rr.start.max(start);
            let r1 = (rr.start + self.nsblk as u64).min(end);
            for s in r0..r1 {
                let src = (s - rr.start) as usize * width;
                let dst = (s - start) as usize * width;
                data[dst..dst + width].copy_from_slice(&vals[src..src + width]);
            }
        }
        Block::new(h, start, n, data)
    }
}
