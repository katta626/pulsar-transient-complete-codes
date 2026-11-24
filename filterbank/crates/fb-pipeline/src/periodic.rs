use fb_defs::{DmPlan, FbError, PeriodicCandidate, Result, RfiMask, Series};
use fb_fft::FftSearchParams;
use fb_ffa::FfaParams;
use fb_fold::{FoldParams, Folded};
use fb_io::Reader;
use rayon::prelude::*;

use crate::single_pulse::RfiConfig;

#[derive(Debug, Clone)]
pub struct PeriodicConfig {
    pub dm_lo: f64,
    pub dm_hi: f64,
    pub dm_step: Option<f64>,
    pub dm_tol: f64,
    pub pulse_width: f64,
    pub dm_batch: usize,
    pub gulp: usize,
    pub fft: Option<FftSearchParams>,
    pub max_accel: f64,
    pub accel_tol: f64,
    pub ffa: Option<FfaParams>,
    pub birdies: Vec<(f64, f64)>,
    pub min_dm_hits: usize,
    pub fold_top: usize,
    pub nbins: usize,
    pub nsub: usize,
    pub rfi: RfiConfig,
    pub external_mask: Option<RfiMask>,
}

impl Default for PeriodicConfig {
    fn default() -> Self {
        Self {
            dm_lo: 0.0,
            dm_hi: 100.0,
            dm_step: None,
            dm_tol: 1.25,
            pulse_width: 1e-3,
            dm_batch: 32,
            gulp: 1 << 18,
            fft: Some(FftSearchParams::default()),
            max_accel: 0.0,
            accel_tol: 1.1,
            ffa: None,
            birdies: vec![],
            min_dm_hits: 1,
            fold_top: 20,
            nbins: 64,
            nsub: 32,
            rfi: RfiConfig::default(),
            external_mask: None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct FoldedCandidate {
    pub candidate: PeriodicCandidate,
    pub profile: Vec<f64>,
    pub reduced_chi2: f64,
    pub folded: Folded,
}

pub struct PeriodicResult {
    pub plan: DmPlan,
    pub mask: RfiMask,
    pub raw: Vec<PeriodicCandidate>,
    pub sifted: Vec<PeriodicCandidate>,
    pub folded: Vec<FoldedCandidate>,
}

pub fn build_mask(reader: &mut dyn Reader, rfi: &RfiConfig, external: Option<&RfiMask>, gulp: usize) -> Result<RfiMask> {
    let h = reader.header();
    let mut mask = match external {
        Some(m) => {
            m.check_applies_to(&h)?;
            m.clone()
        }
        None => RfiMask::empty(&h, rfi.interval, "filterbank"),
    };
    let mut s = 0u64;
    while s < h.nsamples {
        let mut b = reader.read(s, gulp)?;
        if b.nsamp == 0 {
            break;
        }
        if rfi.zero_dm {
            fb_rfi::zero_dm(&mut b, true)?;
        }
        if rfi.interval_stats {
            fb_rfi::interval_stats_mask(&b, rfi.interval_threshold, &mut mask)?;
        }
        if rfi.sk {
            fb_rfi::spectral_kurtosis(&b, rfi.sk_m, 1.0, rfi.sk_sigma, &mut mask)?;
        }
        if rfi.iqrm {
            let r = rfi.iqrm_radius.unwrap_or((b.nchan() / 10).max(1));
            fb_rfi::iqrm_mask(&b, fb_rfi::Statistic::Std, r, rfi.iqrm_threshold, &mut mask)?;
        }
        s += b.nsamp as u64;
    }
    Ok(mask)
}

pub fn dedisperse_all(reader: &mut dyn Reader, plan: &DmPlan, dm_idx: &[usize], mask: &RfiMask, zero_dm: bool, gulp: usize) -> Result<Vec<Series>> {
    let h = reader.header();
    let sub = DmPlan::from_list(plan.freqs.clone(), plan.tsamp, dm_idx.iter().map(|&i| plan.dms[i]).collect())?;
    let maxd = plan.max_delay_samples();
    let total_out = (h.nsamples as usize).saturating_sub(maxd);
    let mut out: Vec<Vec<f32>> = vec![Vec::with_capacity(total_out); sub.ndm()];
    let gulp = gulp.max(2 * maxd).max(4096);
    let kill: Vec<bool> = mask.channel_fraction().iter().map(|f| *f > 0.5).collect();
    let kill = if kill.len() == h.nchan() { Some(kill) } else { None };
    let mut s = 0u64;
    while (s as usize) < total_out {
        let mut b = reader.read(s, gulp + maxd)?;
        if b.nsamp <= maxd {
            break;
        }
        if zero_dm {
            fb_rfi::zero_dm(&mut b, true)?;
        }
        let means: Vec<f32> = fb_rfi::stats::channel_stats(&b.data, b.nsamp, b.nchan()).mean.iter().map(|m| *m as f32).collect();
        b.apply_mask(mask, Some(&means))?;
        let r = fb_dedisp::brute_force(&b, &sub, kill.as_deref())?;
        let take = r.nout.min(gulp).min(total_out - s as usize);
        for d in 0..sub.ndm() {
            out[d].extend_from_slice(&r.trial(d)[..take]);
        }
        s += take as u64;
        if take == 0 {
            break;
        }
    }
    Ok(out
        .into_iter()
        .enumerate()
        .map(|(i, data)| {
            let mut hh = (*h).clone();
            hh.nsamples = data.len() as u64;
            hh.provenance.push(format!("dedisp dm={:.4}", sub.dms[i]));
            Series { header: hh, dm: sub.dms[i], f_ref: plan.f_ref, data }
        })
        .collect())
}

pub fn search_one(s: &Series, cfg: &PeriodicConfig) -> Vec<PeriodicCandidate> {
    let mut c = Vec::new();
    if let Some(p) = &cfg.fft {
        let t_obs = s.len() as f64 * s.tsamp();
        let trials = fb_fft::accel_trials(t_obs, s.tsamp(), cfg.max_accel, cfg.pulse_width, cfg.accel_tol);
        c.extend(fb_fft::search_series(s, &trials, &cfg.birdies, p));
    }
    if let Some(p) = &cfg.ffa {
        c.extend(fb_ffa::search(s, p));
    }
    c
}

pub fn fold_candidate(reader: &mut dyn Reader, c: &PeriodicCandidate, mask: &RfiMask, nbins: usize, nsub: usize, nchan_out: usize, gulp: usize) -> Result<FoldedCandidate> {
    let h = reader.header();
    let fp = FoldParams { f: c.f_hz, fdot: c.fdot, fddot: c.fddot, epoch: c.epoch, dm: c.dm, nbins };
    let mut acc: Option<Folded> = None;
    let mut var_sum = 0f64;
    let mut var_n = 0f64;
    let mut s = 0u64;
    while s < h.nsamples {
        let mut b = reader.read(s, gulp)?;
        if b.nsamp == 0 {
            break;
        }
        let st = fb_rfi::stats::channel_stats(&b.data, b.nsamp, b.nchan());
        let means: Vec<f32> = st.mean.iter().map(|m| *m as f32).collect();
        b.apply_mask(mask, Some(&means))?;
        var_sum += st.std.iter().map(|x| x * x).sum::<f64>();
        var_n += st.std.len() as f64;
        acc = Some(fb_fold::fold_block(&b, &fp, nsub, nchan_out, acc)?);
        s += b.nsamp as u64;
    }
    let folded = acc.ok_or_else(|| FbError::Invalid("no data to fold".into()))?;
    let var = var_sum / var_n.max(1.0);
    let chi = folded.reduced_chi2(var.max(1e-12));
    Ok(FoldedCandidate { candidate: c.clone(), profile: fb_fold::normalised_profile(&folded), reduced_chi2: chi, folded })
}

pub fn run_periodic(reader: &mut dyn Reader, cfg: &PeriodicConfig) -> Result<PeriodicResult> {
    let h = reader.header();
    if h.npol != 1 {
        return Err(FbError::Unsupported("periodicity search needs total intensity".into()));
    }
    let plan = match cfg.dm_step {
        Some(st) => DmPlan::linear(h.freqs.clone(), h.tsamp, cfg.dm_lo, cfg.dm_hi, st)?,
        None => DmPlan::tolerance(h.freqs.clone(), h.tsamp, cfg.dm_lo, cfg.dm_hi, cfg.pulse_width, cfg.dm_tol)?,
    };
    let mask = build_mask(reader, &cfg.rfi, cfg.external_mask.as_ref(), cfg.gulp)?;
    let mut raw = Vec::new();
    let idx: Vec<usize> = (0..plan.ndm()).collect();
    for batch in idx.chunks(cfg.dm_batch.max(1)) {
        let series = dedisperse_all(reader, &plan, batch, &mask, cfg.rfi.zero_dm, cfg.gulp)?;
        let found: Vec<PeriodicCandidate> = series.par_iter().flat_map(|s| search_one(s, cfg)).collect();
        raw.extend(found);
    }
    let t_obs = h.duration();
    let sifted = fb_fft::sift(&raw, t_obs, cfg.min_dm_hits, 16);
    let mut folded = Vec::new();
    for c in sifted.iter().take(cfg.fold_top) {
        folded.push(fold_candidate(reader, c, &mask, cfg.nbins, cfg.nsub, 16, cfg.gulp)?);
    }
    folded.sort_by(|a, b| b.reduced_chi2.partial_cmp(&a.reduced_chi2).unwrap());
    Ok(PeriodicResult { plan, mask, raw, sifted, folded })
}
