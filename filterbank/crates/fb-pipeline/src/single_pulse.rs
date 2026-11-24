use fb_defs::{Block, DmPlan, FbError, Result, RfiMask, SinglePulseCandidate};
use fb_io::Reader;
use fb_rfi::Statistic;
use fb_sp::SearchParams;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Dedisperser {
    Auto,
    BruteForce,
    Subband(usize),
    Fdmt,
}

impl Dedisperser {
    pub fn parse(s: &str, nsub: usize) -> Option<Self> {
        match s {
            "auto" => Some(Dedisperser::Auto),
            "brute" => Some(Dedisperser::BruteForce),
            "subband" => Some(Dedisperser::Subband(nsub)),
            "fdmt" => Some(Dedisperser::Fdmt),
            _ => None,
        }
    }

    pub fn resolve(self, nchan: usize) -> Self {
        match self {
            Dedisperser::Auto if nchan >= 256 => Dedisperser::Subband(((nchan as f64).sqrt() as usize).next_power_of_two().max(16)),
            Dedisperser::Auto => Dedisperser::BruteForce,
            d => d,
        }
    }

    pub fn run(self, block: &Block, plan: &DmPlan, kill: Option<&[bool]>) -> Result<fb_dedisp::DedispOutput> {
        match self.resolve(block.nchan()) {
            Dedisperser::Subband(n) => fb_dedisp::subband(block, plan, n, kill),
            Dedisperser::Fdmt => fb_dedisp::fdmt(block, plan, kill),
            _ => fb_dedisp::brute_force(block, plan, kill),
        }
    }
}

#[derive(Debug, Clone)]
pub struct RfiConfig {
    pub zero_dm: bool,
    pub iqrm: bool,
    pub iqrm_radius: Option<usize>,
    pub iqrm_threshold: f64,
    pub sk: bool,
    pub sk_m: usize,
    pub sk_sigma: f64,
    pub interval_stats: bool,
    pub interval_threshold: f64,
    pub interval: f64,
}

impl Default for RfiConfig {
    fn default() -> Self {
        Self {
            zero_dm: false,
            iqrm: true,
            iqrm_radius: None,
            iqrm_threshold: 3.0,
            sk: false,
            sk_m: 256,
            sk_sigma: 5.0,
            interval_stats: false,
            interval_threshold: 6.0,
            interval: 1.0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct SinglePulseConfig {
    pub dm_lo: f64,
    pub dm_hi: f64,
    pub dm_step: Option<f64>,
    pub dm_tol: f64,
    pub pulse_width: f64,
    pub gulp: usize,
    pub search: SearchParams,
    pub rfi: RfiConfig,
    pub dedisperser: Dedisperser,
    pub external_mask: Option<RfiMask>,
}

impl Default for SinglePulseConfig {
    fn default() -> Self {
        Self {
            dm_lo: 0.0,
            dm_hi: 1000.0,
            dm_step: None,
            dm_tol: 1.25,
            pulse_width: 1e-3,
            gulp: 1 << 16,
            search: SearchParams::default(),
            rfi: RfiConfig::default(),
            dedisperser: Dedisperser::Auto,
            external_mask: None,
        }
    }
}

pub struct SinglePulseResult {
    pub plan: DmPlan,
    pub mask: RfiMask,
    pub candidates: Vec<SinglePulseCandidate>,
}

pub fn make_plan(reader: &dyn Reader, cfg: &SinglePulseConfig) -> Result<DmPlan> {
    let h = reader.header();
    match cfg.dm_step {
        Some(step) => DmPlan::linear(h.freqs.clone(), h.tsamp, cfg.dm_lo, cfg.dm_hi, step),
        None => DmPlan::tolerance(h.freqs.clone(), h.tsamp, cfg.dm_lo, cfg.dm_hi, cfg.pulse_width, cfg.dm_tol),
    }
}

fn clean_block(block: &mut Block, cfg: &RfiConfig, mask: &mut RfiMask) -> Result<Vec<bool>> {
    if cfg.zero_dm {
        fb_rfi::zero_dm(block, true)?;
    }
    if cfg.interval_stats {
        fb_rfi::interval_stats_mask(block, cfg.interval_threshold, mask)?;
    }
    if cfg.sk {
        fb_rfi::spectral_kurtosis(block, cfg.sk_m, 1.0, cfg.sk_sigma, mask)?;
    }
    if cfg.iqrm {
        let r = cfg.iqrm_radius.unwrap_or((block.nchan() / 10).max(1));
        fb_rfi::iqrm_mask(block, Statistic::Std, r, cfg.iqrm_threshold, mask)?;
    }
    let segs = mask.to_indices(&block.header, block.start_sample, block.nsamp)?;
    let mut kill = vec![false; block.nchan()];
    let nchan = block.nchan();
    for s in &segs {
        if (s.t1 - s.t0) * 2 >= block.nsamp {
            for c in 0..nchan {
                kill[c] |= s.bad[c];
            }
        }
    }
    let means = fb_rfi::stats::channel_stats(&block.data, block.nsamp, nchan).mean;
    let fill: Vec<f32> = means.iter().map(|m| *m as f32).collect();
    block.apply_mask(mask, Some(&fill))?;
    Ok(kill)
}

pub fn run_single_pulse(reader: &mut dyn Reader, cfg: &SinglePulseConfig) -> Result<SinglePulseResult> {
    let h = reader.header();
    if h.npol != 1 {
        return Err(FbError::Unsupported("single-pulse search needs total intensity; open PSRFITS with sum_pols".into()));
    }
    let plan = make_plan(reader, cfg)?;
    let maxd = plan.max_delay_samples();
    let gulp = cfg.gulp.max(2 * maxd).max(4096);
    let mut mask = match &cfg.external_mask {
        Some(m) => {
            m.check_applies_to(&h)?;
            m.clone()
        }
        None => RfiMask::empty(&h, cfg.rfi.interval, "filterbank"),
    };
    let ws = fb_sp::widths(cfg.search.max_width, cfg.search.snr_loss);
    let mut all = Vec::new();
    let mut start = 0u64;
    while start < h.nsamples {
        let mut block = reader.read(start, gulp + maxd)?;
        if block.nsamp <= maxd + 16 {
            break;
        }
        let kill = clean_block(&mut block, &cfg.rfi, &mut mask)?;
        let out = cfg.dedisperser.run(&block, &plan, Some(&kill))?;
        let trials: Vec<Vec<f32>> = (0..out.ndm).map(|i| out.trial(i).to_vec()).collect();
        let events = fb_sp::search_trials(&trials, &cfg.search);
        let clusters = fb_sp::cluster(&events, &cfg.search, &ws);
        let cands = fb_sp::to_candidates(&clusters, &h, &block.start_epoch(), block.start_sample, &plan.dms, plan.f_ref);
        all.extend(cands);
        start += gulp as u64;
    }
    all.sort_by(|a, b| a.sample.cmp(&b.sample));
    Ok(SinglePulseResult { plan, mask, candidates: all })
}
