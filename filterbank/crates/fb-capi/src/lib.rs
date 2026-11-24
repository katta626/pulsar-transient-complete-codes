use std::cell::RefCell;
use std::ffi::{c_char, CStr, CString};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use fb_defs::{Block, DmPlan, FbError, ObsHeader, Result, RfiMask, Series};
use fb_formats::presto::{read_dat, write_dat, PrestoMask};
use fb_io::{open_with, PsrfitsOptions, Reader, ScalePolicy, SigprocWriter, WriteOptions, ZeroOffPolicy};
use fb_pipeline::{run_single_pulse, Dedisperser, SinglePulseConfig};
use serde_json::{json, Value};

thread_local! {
    static LAST_ERROR: RefCell<CString> = RefCell::new(CString::new("").unwrap());
}

fn set_error(e: impl std::fmt::Display) {
    LAST_ERROR.with(|c| *c.borrow_mut() = CString::new(e.to_string().replace('\0', " ")).unwrap());
}

fn cstr<'a>(p: *const c_char) -> Result<&'a str> {
    if p.is_null() {
        return Err(FbError::Invalid("null string".into()));
    }
    unsafe { CStr::from_ptr(p) }.to_str().map_err(|e| FbError::Invalid(e.to_string()))
}

fn out_string(s: String) -> *mut c_char {
    CString::new(s).map(|c| c.into_raw()).unwrap_or(std::ptr::null_mut())
}

fn jparse<T: serde::de::DeserializeOwned>(p: *const c_char) -> Result<T> {
    Ok(serde_json::from_str(cstr(p)?)?)
}

fn jstr<T: serde::Serialize>(v: &T) -> Result<*mut c_char> {
    Ok(out_string(serde_json::to_string(v)?))
}

macro_rules! guard_ptr {
    ($e:expr) => {
        match (|| -> Result<*mut c_char> { $e })() {
            Ok(p) => p,
            Err(e) => {
                set_error(e);
                std::ptr::null_mut()
            }
        }
    };
}

macro_rules! guard_i64 {
    ($e:expr) => {
        match (|| -> Result<i64> { $e })() {
            Ok(v) => v,
            Err(e) => {
                set_error(e);
                -1
            }
        }
    };
}

#[no_mangle]
pub extern "C" fn fb_version() -> *const c_char {
    concat!(env!("CARGO_PKG_VERSION"), "\0").as_ptr() as *const c_char
}

#[no_mangle]
pub extern "C" fn fb_last_error() -> *const c_char {
    LAST_ERROR.with(|c| c.borrow().as_ptr())
}

#[no_mangle]
pub extern "C" fn fb_free_string(p: *mut c_char) {
    if !p.is_null() {
        unsafe { drop(CString::from_raw(p)) };
    }
}

pub struct Handle {
    reader: Box<dyn Reader>,
}

fn psrfits_opts(v: &Value) -> PsrfitsOptions {
    let mut o = PsrfitsOptions::default();
    if let Some(z) = v.get("zero_off").and_then(|x| x.as_str()) {
        o.zero_off = match z {
            "dspsr" => ZeroOffPolicy::Dspsr,
            "ignore" => ZeroOffPolicy::Ignore,
            _ => ZeroOffPolicy::Spec,
        };
    }
    if let Some(s) = v.get("scales").and_then(|x| x.as_str()) {
        o.scales = match s {
            "always" => ScalePolicy::Always,
            "never" => ScalePolicy::Never,
            _ => ScalePolicy::Auto,
        };
    }
    if let Some(b) = v.get("weights").and_then(|x| x.as_bool()) {
        o.weights = b;
    }
    if let Some(b) = v.get("sum_pols").and_then(|x| x.as_bool()) {
        o.sum_pols = b;
    }
    if v.get("raw").and_then(|x| x.as_bool()) == Some(true) {
        o = PsrfitsOptions::raw();
    }
    o
}

#[no_mangle]
pub extern "C" fn fb_open(paths_json: *const c_char, opts_json: *const c_char) -> *mut Handle {
    let r = (|| -> Result<Box<Handle>> {
        let paths: Vec<PathBuf> = jparse(paths_json)?;
        let opts: Value = if opts_json.is_null() { json!({}) } else { jparse(opts_json)? };
        Ok(Box::new(Handle { reader: open_with(&paths, psrfits_opts(&opts))? }))
    })();
    match r {
        Ok(h) => Box::into_raw(h),
        Err(e) => {
            set_error(e);
            std::ptr::null_mut()
        }
    }
}

#[no_mangle]
pub extern "C" fn fb_close(h: *mut Handle) {
    if !h.is_null() {
        unsafe { drop(Box::from_raw(h)) };
    }
}

#[no_mangle]
pub extern "C" fn fb_header(h: *mut Handle) -> *mut c_char {
    guard_ptr!({
        let h = unsafe { h.as_ref() }.ok_or_else(|| FbError::Invalid("null handle".into()))?;
        let gaps: Vec<(u64, u64)> = h.reader.gaps().iter().map(|g| (g.start, g.len)).collect();
        let mut v = serde_json::to_value(&*h.reader.header())?;
        v["gaps"] = json!(gaps);
        jstr(&v)
    })
}

#[no_mangle]
pub extern "C" fn fb_read(h: *mut Handle, start: u64, nsamp: u64, out: *mut f32, cap: u64) -> i64 {
    guard_i64!({
        let h = unsafe { h.as_mut() }.ok_or_else(|| FbError::Invalid("null handle".into()))?;
        let b = h.reader.read(start, nsamp as usize)?;
        if b.data.len() as u64 > cap {
            return Err(FbError::Invalid(format!("buffer holds {cap} values, need {}", b.data.len())));
        }
        unsafe { std::ptr::copy_nonoverlapping(b.data.as_ptr(), out, b.data.len()) };
        Ok(b.nsamp as i64)
    })
}

fn block_from(header_json: *const c_char, data: *const f32, nsamp: u64, start_sample: u64) -> Result<Block> {
    let h: ObsHeader = jparse(header_json)?;
    let n = nsamp as usize * h.npol * h.nchan();
    let v = unsafe { std::slice::from_raw_parts(data, n) }.to_vec();
    Block::new(Arc::new(h), start_sample, nsamp as usize, v)
}

#[no_mangle]
pub extern "C" fn fb_write_fil(path: *const c_char, header_json: *const c_char, nbits: u32, data: *const f32, nvalues: u64) -> i64 {
    guard_i64!({
        let h: ObsHeader = jparse(header_json)?;
        let p = PathBuf::from(cstr(path)?);
        let mut w = SigprocWriter::create(&p, &h, WriteOptions { nbits, data_type: 1, refdm: None, rawdatafile: h.source_name.clone() })?;
        w.write_values(unsafe { std::slice::from_raw_parts(data, nvalues as usize) })?;
        Ok(w.finish()? as i64)
    })
}

#[no_mangle]
pub extern "C" fn fb_plan(header_json: *const c_char, lo: f64, hi: f64, step: f64, width: f64, tol: f64) -> *mut c_char {
    guard_ptr!({
        let h: ObsHeader = jparse(header_json)?;
        let p = if step > 0.0 {
            DmPlan::linear(h.freqs.clone(), h.tsamp, lo, hi, step)?
        } else {
            DmPlan::tolerance(h.freqs.clone(), h.tsamp, lo, hi, width, tol)?
        };
        jstr(&p)
    })
}

#[no_mangle]
pub extern "C" fn fb_plan_from_list(header_json: *const c_char, dms_json: *const c_char) -> *mut c_char {
    guard_ptr!({
        let h: ObsHeader = jparse(header_json)?;
        let dms: Vec<f64> = jparse(dms_json)?;
        jstr(&DmPlan::from_list(h.freqs.clone(), h.tsamp, dms)?)
    })
}

#[no_mangle]
pub extern "C" fn fb_delay_table(plan_json: *const c_char, out: *mut i32, cap: u64) -> i64 {
    guard_i64!({
        let p: DmPlan = jparse(plan_json)?;
        let t = p.delay_table();
        if t.len() as u64 > cap {
            return Err(FbError::Invalid("delay buffer too small".into()));
        }
        unsafe { std::ptr::copy_nonoverlapping(t.as_ptr(), out, t.len()) };
        Ok(p.max_delay_samples() as i64)
    })
}

#[no_mangle]
pub extern "C" fn fb_dedisperse(
    header_json: *const c_char,
    plan_json: *const c_char,
    data: *const f32,
    nsamp: u64,
    killmask: *const u8,
    nsub: u32,
    out: *mut f32,
    cap: u64,
) -> i64 {
    guard_i64!({
        let b = block_from(header_json, data, nsamp, 0)?;
        let p: DmPlan = jparse(plan_json)?;
        let kill: Option<Vec<bool>> = if killmask.is_null() {
            None
        } else {
            Some(unsafe { std::slice::from_raw_parts(killmask, b.nchan()) }.iter().map(|v| *v != 0).collect())
        };
        let r = if nsub == u32::MAX {
            fb_dedisp::fdmt(&b, &p, kill.as_deref())?
        } else if nsub > 0 {
            fb_dedisp::subband(&b, &p, nsub as usize, kill.as_deref())?
        } else {
            fb_dedisp::brute_force(&b, &p, kill.as_deref())?
        };
        if r.data.len() as u64 > cap {
            return Err(FbError::Invalid(format!("output holds {cap} values, need {}", r.data.len())));
        }
        unsafe { std::ptr::copy_nonoverlapping(r.data.as_ptr(), out, r.data.len()) };
        Ok(r.nout as i64)
    })
}

#[no_mangle]
pub extern "C" fn fb_iqrm(x: *const f64, n: u64, radius: u32, threshold: f64, out: *mut u8) -> i64 {
    guard_i64!({
        let v = unsafe { std::slice::from_raw_parts(x, n as usize) };
        let m = fb_rfi::iqrm(v, radius as usize, threshold);
        let o = unsafe { std::slice::from_raw_parts_mut(out, n as usize) };
        for (a, b) in o.iter_mut().zip(&m) {
            *a = *b as u8;
        }
        Ok(m.iter().filter(|b| **b).count() as i64)
    })
}

#[no_mangle]
pub extern "C" fn fb_mask_new(header_json: *const c_char, interval: f64) -> *mut c_char {
    guard_ptr!({
        let h: ObsHeader = jparse(header_json)?;
        jstr(&RfiMask::empty(&h, interval, "filterbank"))
    })
}

#[no_mangle]
pub extern "C" fn fb_rfi_block(
    header_json: *const c_char,
    data: *mut f32,
    nsamp: u64,
    start_sample: u64,
    mask_json: *const c_char,
    cfg_json: *const c_char,
) -> *mut c_char {
    guard_ptr!({
        let mut b = block_from(header_json, data, nsamp, start_sample)?;
        let mut m: RfiMask = jparse(mask_json)?;
        let c: Value = if cfg_json.is_null() { json!({}) } else { jparse(cfg_json)? };
        let f = |k: &str, d: f64| c.get(k).and_then(|v| v.as_f64()).unwrap_or(d);
        let on = |k: &str, d: bool| c.get(k).and_then(|v| v.as_bool()).unwrap_or(d);
        if on("zero_dm", false) {
            fb_rfi::zero_dm(&mut b, true)?;
        }
        if on("stats", false) {
            fb_rfi::interval_stats_mask(&b, f("stats_threshold", 6.0), &mut m)?;
        }
        if on("sk", false) {
            fb_rfi::spectral_kurtosis(&b, f("sk_m", 256.0) as usize, f("sk_n", 1.0), f("sk_sigma", 5.0), &mut m)?;
        }
        if on("iqrm", true) {
            let r = c.get("iqrm_radius").and_then(|v| v.as_u64()).map(|v| v as usize).unwrap_or((b.nchan() / 10).max(1));
            fb_rfi::iqrm_mask(&b, fb_rfi::Statistic::Std, r, f("iqrm_threshold", 3.0), &mut m)?;
        }
        if on("zero_dm", false) {
            unsafe { std::ptr::copy_nonoverlapping(b.data.as_ptr(), data, b.data.len()) };
        }
        jstr(&m)
    })
}

#[no_mangle]
pub extern "C" fn fb_mask_apply(header_json: *const c_char, data: *mut f32, nsamp: u64, start_sample: u64, mask_json: *const c_char, fill_mean: u8) -> i64 {
    guard_i64!({
        let mut b = block_from(header_json, data, nsamp, start_sample)?;
        let m: RfiMask = jparse(mask_json)?;
        let fill: Option<Vec<f32>> = if fill_mean != 0 {
            Some(fb_rfi::stats::channel_stats(&b.data, b.nsamp, b.nchan() * b.npol()).mean.iter().map(|v| *v as f32).collect())
        } else {
            None
        };
        let n = b.apply_mask(&m, fill.as_deref())?;
        unsafe { std::ptr::copy_nonoverlapping(b.data.as_ptr(), data, b.data.len()) };
        Ok(n as i64)
    })
}

#[no_mangle]
pub extern "C" fn fb_presto_mask_read(path: *const c_char, header_json: *const c_char) -> *mut c_char {
    guard_ptr!({
        let pm = PrestoMask::read(Path::new(cstr(path)?))?;
        let m = if header_json.is_null() {
            pm.to_mask("unknown")
        } else {
            pm.to_mask_for(&jparse::<ObsHeader>(header_json)?)
        };
        jstr(&m)
    })
}

#[no_mangle]
pub extern "C" fn fb_presto_mask_write(mask_json: *const c_char, tsamp: f64, path: *const c_char) -> i64 {
    guard_i64!({
        let m: RfiMask = jparse(mask_json)?;
        PrestoMask::from_mask(&m, tsamp)?.write(Path::new(cstr(path)?))?;
        Ok(0)
    })
}

#[no_mangle]
pub extern "C" fn fb_dat_info(path: *const c_char) -> *mut c_char {
    guard_ptr!({
        let s = read_dat(Path::new(cstr(path)?))?;
        jstr(&json!({"header": s.header, "dm": s.dm, "f_ref": s.f_ref, "n": s.data.len()}))
    })
}

#[no_mangle]
pub extern "C" fn fb_dat_read(path: *const c_char, out: *mut f32, cap: u64) -> i64 {
    guard_i64!({
        let s = read_dat(Path::new(cstr(path)?))?;
        if s.data.len() as u64 > cap {
            return Err(FbError::Invalid("buffer too small".into()));
        }
        unsafe { std::ptr::copy_nonoverlapping(s.data.as_ptr(), out, s.data.len()) };
        Ok(s.data.len() as i64)
    })
}

#[no_mangle]
pub extern "C" fn fb_dat_write(root: *const c_char, header_json: *const c_char, dm: f64, f_ref: f64, data: *const f32, n: u64) -> i64 {
    guard_i64!({
        let header: ObsHeader = jparse(header_json)?;
        let s = Series { header, dm, f_ref, data: unsafe { std::slice::from_raw_parts(data, n as usize) }.to_vec() };
        write_dat(Path::new(cstr(root)?), &s)?;
        Ok(n as i64)
    })
}

#[no_mangle]
pub extern "C" fn fb_sp_search(
    header_json: *const c_char,
    plan_json: *const c_char,
    trials: *const f32,
    ndm: u64,
    nout: u64,
    start_sample: u64,
    params_json: *const c_char,
) -> *mut c_char {
    guard_ptr!({
        let h: ObsHeader = jparse(header_json)?;
        let p: DmPlan = jparse(plan_json)?;
        let c: Value = if params_json.is_null() { json!({}) } else { jparse(params_json)? };
        let mut sp = fb_sp::SearchParams::default();
        if let Some(v) = c.get("threshold").and_then(|v| v.as_f64()) {
            sp.threshold = v as f32;
        }
        if let Some(v) = c.get("max_width").and_then(|v| v.as_u64()) {
            sp.max_width = v as usize;
        }
        let all = unsafe { std::slice::from_raw_parts(trials, (ndm * nout) as usize) };
        let tr: Vec<Vec<f32>> = all.chunks(nout as usize).map(|x| x.to_vec()).collect();
        let ws = fb_sp::widths(sp.max_width, sp.snr_loss);
        let ev = fb_sp::search_trials(&tr, &sp);
        let cl = fb_sp::cluster(&ev, &sp, &ws);
        let start = h.epoch_of_sample(start_sample);
        jstr(&fb_sp::to_candidates(&cl, &h, &start, start_sample, &p.dms, p.f_ref))
    })
}

#[no_mangle]
pub extern "C" fn fb_run_sp(paths_json: *const c_char, cfg_json: *const c_char) -> *mut c_char {
    guard_ptr!({
        let paths: Vec<PathBuf> = jparse(paths_json)?;
        let c: Value = if cfg_json.is_null() { json!({}) } else { jparse(cfg_json)? };
        let mut reader = open_with(&paths, psrfits_opts(&c))?;
        let f = |k: &str| c.get(k).and_then(|v| v.as_f64());
        let b = |k: &str| c.get(k).and_then(|v| v.as_bool());
        let mut cfg = SinglePulseConfig::default();
        if let Some(v) = f("dm_lo") {
            cfg.dm_lo = v;
        }
        if let Some(v) = f("dm_hi") {
            cfg.dm_hi = v;
        }
        cfg.dm_step = f("dm_step");
        if let Some(v) = f("dm_tol") {
            cfg.dm_tol = v;
        }
        if let Some(v) = f("width_s") {
            cfg.pulse_width = v;
        }
        if let Some(v) = f("threshold") {
            cfg.search.threshold = v as f32;
        }
        if let Some(v) = f("max_width") {
            cfg.search.max_width = v as usize;
        }
        if let Some(v) = f("gulp") {
            cfg.gulp = v as usize;
        }
        if let Some(v) = f("subband") {
            cfg.dedisperser = Dedisperser::Subband(v as usize);
        }
        if let Some(m) = c.get("dedisp").and_then(|v| v.as_str()) {
            cfg.dedisperser = Dedisperser::parse(m, f("subband").unwrap_or(64.0) as usize).ok_or_else(|| FbError::Invalid(format!("dedisp '{m}'")))?;
        }
        if let Some(v) = b("zero_dm") {
            cfg.rfi.zero_dm = v;
        }
        if let Some(v) = b("iqrm") {
            cfg.rfi.iqrm = v;
        }
        if let Some(v) = b("sk") {
            cfg.rfi.sk = v;
        }
        if let Some(v) = b("stats") {
            cfg.rfi.interval_stats = v;
        }
        if let Some(v) = f("interval") {
            cfg.rfi.interval = v;
        }
        if let Some(m) = c.get("mask") {
            if !m.is_null() {
                cfg.external_mask = Some(serde_json::from_value(m.clone())?);
            }
        }
        let r = run_single_pulse(reader.as_mut(), &cfg)?;
        jstr(&json!({"plan": r.plan, "mask": r.mask, "candidates": r.candidates}))
    })
}
