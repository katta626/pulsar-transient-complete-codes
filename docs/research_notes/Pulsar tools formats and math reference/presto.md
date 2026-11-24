# PRESTO (Scott Ransom's pulsar search package): programs, workflows, file formats, algorithms, versions

Source basis: a shallow clone of https://github.com/scottransom/presto (default branch `master`, HEAD commit `4a5bd80`, 2026-11-18, "Fix prepfold on short .dat files and pgopen() on aarch64"; i.e. development after v6.3.3), read directly: `src/*.c`, `include/*.h`, `clig/*_cmd.cli` (the option specs, which hold the defaults), `bin/*.py`, `python/presto/*.py`, `docs/*.1` man pages, `docs/PRESTO_search_tutorial.pdf`, `CHANGELOG.md`, `FAQ.md`, `ROADMAP.md`, `COPYING`, `CITATION.cff`. Release dates come from the GitHub releases API. Paper abstracts come from arXiv.
Short link convention: `B/` = `https://github.com/scottransom/presto/blob/master/`.

---

## 1. Catalogue of the main PRESTO programs: inputs, outputs, header/data fields used, algorithm, source files

### Takeaway
PRESTO is a GPL-2.0 C library (`libpresto`) with C command-line tools and Python scripts. Raw data (SIGPROC filterbank or PSRFITS) go through `rfifind` (mask), `prepdata`/`prepsubband`/`mpiprepsubband` (dedispersed, usually barycentred `.dat` + `.inf`), `realfft` (`.fft`), `rednoise` (`_red.fft`), `zapbirds`/`simple_zapbirds.py` (birdies), and then `accelsearch` (`_ACCEL_<z>[_JERK_<w>]` + `.cand`) or, new in v6.2, `coherentsearch` (`.cohout`). After that come `ACCEL_sift.py` and `prepfold` (`.pfd`, `.bestprof`, `.ps`). Single pulses go `.dat` → `single_pulse_search.py` → `.singlepulse`, and timing goes `.pfd` → `get_TOAs.py`. There is no FFA program in PRESTO. `search_bin` (phase-modulation/sideband search) exists but is rarely used.

### Cited Findings

**readfile**
- Prints the header of a raw data file (filterbank or PSRFITS), e.g. "Spectra per subint". The tutorial's step 1 is "Examine data format (readfile)". — [PRESTO_search_tutorial.pdf](B/docs/PRESTO_search_tutorial.pdf); source `src/readfile.c`, `clig/readfile_cmd.cli`. For SIGPROC files the block size ("spectra per subint") is hard-coded as `s->spectra_per_subint = 2400` in `src/sigproc_fb.c`. For PSRFITS it is the SUBINT row length. `-downsamp` must divide it. — [FAQ.md](B/FAQ.md) ("How does PRESTO decide what block-size...", "...not letting me use certain downsampling factors").

**rfifind** (`src/rfifind.c`, `src/rfifind_plot.c`, `src/rfi_utils.c`, `src/mask.c`, `src/minifft.c`, `clig/rfifind_cmd.cli`, `include/rfifind.h`)
- Purpose: it "divides the observation into ... second (or spectra) chunks, computes time-domain and Fourier-domain statistics for each chunk of each channel, and flags outliers". A whole interval is masked if more than `-chanfrac` of its channels are bad, and a whole channel if more than `-intfrac` of its intervals are bad. — [docs/rfifind.1](B/docs/rfifind.1)
- Defaults (clig): `-time 30.0` s, `-timesig 10`, `-freqsig 4`, `-chanfrac 0.7`, `-intfrac 0.3`, `-clip 6.0` (time-domain clipping always on unless `-noclip`). Also `-zerodm`, `-zapchan`, `-zapints`, `-ignorechan`, `-mask` (re-use a mask), `-nocompute` (re-derive the mask from the existing `.stats` with new thresholds, without reading raw data), `-ncpus` (OpenMP, v5.3). — [clig/rfifind_cmd.cli](B/clig/rfifind_cmd.cli), [docs/rfifind.1](B/docs/rfifind.1)
- Outputs: `<o>_rfifind.mask`, `.stats`, `.bytemask`, `.rfi`, plus `.inf` and `.ps` plots. — [src/rfifind.c](B/src/rfifind.c) (the `sprintf` calls for `_rfifind`, `.mask`, `.rfi`, `.stats`, `.bytemask`)
- Per (interval ii, channel jj) it stores `dataavg` (mean), `datastd` (sqrt of variance via `avg_var`) and `datapow` (the max normalized power from a real FFT of that chunk, searched with `search_fft()` using `RFI_NUMHARMSUM 4`, `RFI_NUMBETWEEN 2` (interbinning), `RFI_LOBIN 5`). The FFT normalization is `norm = std^2 * ptsperint`. Chunks with std ≤ 1e-4 are not FFT-searched. Padding chunks get avg 0, std 0, pow 1. Detected birdies are merged into an RFI list when the fractional frequency match is within `RFI_FRACTERROR 0.002`. — [src/rfifind.c](B/src/rfifind.c), [include/rfifind.h](B/include/rfifind.h)
- `.stats` binary layout: header ints/floats, then `datapow`, `dataavg`, `datastd` as float32 arrays of `numchan*numint`. — [src/rfifind.c](B/src/rfifind.c) `write_statsfile`/`read_statsfile`
- `.mask` binary layout: doubles `timesigma, freqsigma, mjd, dtint, lofreq, dfreq`; ints `numchan, numint, ptsperint, num_zap_chans`, then `zap_chans[]`, `num_zap_ints`, `zap_ints[]`, `num_chans_per_int[numint]`, then per-interval channel lists. — [src/mask.c](B/src/mask.c) `read_mask`/`write_mask`
- FAQ: rfifind cannot remove mains 50/60 Hz from accelsearch results (that needs birdie zapping), and it uses the same channel numbering as the raw file. — [FAQ.md](B/FAQ.md) RFI section

**DDplan.py** (`bin/DDplan.py`)
- Plans DM steps, downsampling and subband DM steps so that total smearing stays under a tolerance. `total_smear = sqrt(2*(1000 dt)^2 + chan_smear(DM)^2 + subband_smear^2 + BW_smear(dDM)^2)` in ms. `BW_smear` uses a DM error of 0.5*dDM over the full BW. `subband_smear` uses 0.5*subDMstep over BW/numsub. `guess_DMstep_approx = dt*0.0001205*f_ctr^3/(0.5*BW)`. Allowed dDMs are `[0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 300]`, with fudge factor `ff = 1.2`. Downsampling grows when channel smearing exceeds dt. — [bin/DDplan.py](B/bin/DDplan.py)
- It can read observation parameters from a filterbank/PSRFITS file, and `-w` writes a `dedisp_*.py` script (since v3.0.1). Since v5.1.0 it handles DM smearing more correctly. — [CHANGELOG.md](B/CHANGELOG.md); tutorial example `DDplan.py -d 500.0 -n 96 -b 96 -t 0.000072 -f 1400.0 -s 32 -r 0.5` — [tutorial](B/docs/PRESTO_search_tutorial.pdf)

**prepdata** (`src/prepdata.c`, `src/dispersion.c`, `src/barycenter.c`, `src/observatories.c`)
- Dedisperses raw data at one DM (or reads a time series) and writes `.dat` (float32) plus `.inf`. By default it barycentres; `-nobary` turns that off. — [tutorial](B/docs/PRESTO_search_tutorial.pdf), [docs/prepdata.1](B/docs/prepdata.1)
- Barycentring: topocentric TOAs are computed every `TDT = 20.0` s and converted with `barycenter()` to barycentric TDB MJDs plus v/c. The average v/c Doppler-corrects the channel frequencies (`doppler(tobsf, avgvoverc)`). The bary−topo differences, in output bins, become a list of `diffbins` where single samples are added (positive) or removed (negative), so the output is on a uniform barycentric grid. — [src/prepdata.c](B/src/prepdata.c) (lines with `TDT`, `numbarypts`, `diffbins`)
- Since v6.0.0 `barycenter()` uses ERFA (built-in analytic ephemeris, leap seconds, internal observatory table) instead of TEMPO. Agreement is about 10 µs absolute and about 1 µs differential, the `ephem` argument (e.g. DE405) is ignored, and the old code is kept as `barycenter_tempo()`. — [CHANGELOG.md](B/CHANGELOG.md) v6.0.0; [src/barycenter.c](B/src/barycenter.c)

**prepsubband / mpiprepsubband** (`src/prepsubband.c`, `src/mpiprepsubband.c`, `src/mpiprepsubband_utils.c`, `src/dispersion.c`)
- Makes many dedispersed `.dat`/`.inf` pairs in one pass. It always uses two stages: channels → `-nsub` subbands partially dedispersed to a nominal DM, then subbands → each trial DM. `-sub -subdm X` writes persistent subband files `<o>_DM<X>.subNN` for later re-dedispersion or folding. — [docs/prepsubband.1](B/docs/prepsubband.1); format string `"%s_DM%.*f.sub%0Nd"` in [src/prepsubband.c](B/src/prepsubband.c)
- Defaults: `-numdms 10`, `-dmstep 1.0`, `-lodm 0`, `-nsub 32`, `-downsamp 1` (max 128), `-dmprec 2` (DM decimals in filenames). Also `-mask`, `-ignorechan`, `-zerodm`, `-nobary`, `-numout`, `-offset`/`-start`. — [clig/prepsubband_cmd.cli](B/clig/prepsubband_cmd.cli)
- mpiprepsubband: one CPU reads the data and broadcasts it to the other CPUs, each of which effectively runs one DDplan "call". — [tutorial](B/docs/PRESTO_search_tutorial.pdf); [FAQ.md](B/FAQ.md) mpiprepsubband entry
- Delays: `delay_from_dm(dm,f) = dm/(0.000241 f^2)` s (f in MHz). Channel frequency is the channel centre. The highest-frequency channel is never shifted; lower channels are brought "forward". The `.inf` epoch is the start of the first sample at the centre of the highest channel. — [src/dispersion.c](B/src/dispersion.c); [FAQ.md](B/FAQ.md) "where in a frequency channel..."
- `subband_search_delays()` computes per-channel delays at the DM and then subtracts the delay of each subband's highest channel. Each subband is therefore internally dedispersed, while the subbands stay offset from one another so that `float_dedisp()` can finish the job. — [src/dispersion.c](B/src/dispersion.c)
- v5.3.0 turned off OpenMP in the dedispersion routines. — [CHANGELOG.md](B/CHANGELOG.md)

**realfft** (`src/realfft.c`, `src/fastffts.c`, `src/fftcalls.c`, `src/twopass*.c`)
- Forward FFT `.dat` → `.fft`, or inverse `.fft` → `.dat`, with the direction set by the suffix (`-fwd`/`-inv` override). It accepts many files per call (v4). Above about 1e9 samples it switches to an out-of-core two-pass FFT (`-tmpdir`). — [clig/realfft_cmd.cli](B/clig/realfft_cmd.cli)
- Packing: N real samples → N/2 complex64 values. Element 0's real part holds DC and its imaginary part holds the Nyquist value. The forward sign is −1. It is a length-N/2 complex FFT followed by the Numerical-Recipes `realft` untangling, with a sign fix for the n/2 imaginary part. — [src/fastffts.c](B/src/fastffts.c)
- `realfft -rednoise` (v6.1.0) writes `<root>_red.fft`/`_red.inf` directly, bit-identical to running realfft then rednoise. — [CHANGELOG.md](B/CHANGELOG.md)

**rednoise** (`src/rednoise.c`, `deredden()`/`deredden_file()` in `src/misc_utils.c`)
- Defaults: `-startwidth 6`, `-endwidth 100`, `-endfreq 6` Hz. `-inplace` (v6.1.0) whitens in place and renames to `_red.fft`. It writes `_red.inf` (v3.0.1+). — [clig/rednoise_cmd.cli](B/clig/rednoise_cmd.cli), [docs/rednoise.1](B/docs/rednoise.1)

**zapbirds / simple_zapbirds.py / makezaplist.py / accelsearch -zaplist** (`src/zapbirds.c`, `src/birdzap.c`, `bin/simple_zapbirds.py`, `bin/makezaplist.py`)
- A `.birds` file has 5 ASCII columns: freq (Hz), width (Hz), number of harmonics, grow (0/1, whether the width grows with harmonic number), bary (0/1). `#` starts a comment. `makezaplist.py X.birds` (it needs `X.inf`) → `X.zaplist`, then `zapbirds -zap -zapfile X.zaplist -baryv <v/c> file.fft`. `simple_zapbirds.py X.birds *.fft` does all of this, computing baryv itself. — [tutorial](B/docs/PRESTO_search_tutorial.pdf); [docs/zapbirds.1](B/docs/zapbirds.1). makezaplist zaps at least 40 bins for known pulsars. — [bin/makezaplist.py](B/bin/makezaplist.py)
- Since v6.1.0, zapped bins are filled with χ²(2 dof) noise scaled to the local mean power, with random phases and a fixed seed, instead of a constant amplitude. `-constamp`/`--constamp` gives the old behaviour. The local level is measured 100 bins clear of the region. — [CHANGELOG.md](B/CHANGELOG.md) v6.1.0; [docs/zapbirds.1](B/docs/zapbirds.1)
- Birdies are found by running `accelsearch -numharm 4 -zmax 0` on the DM 0 topocentric series and inspecting with `explorefft`. — [tutorial](B/docs/PRESTO_search_tutorial.pdf)

**accelsearch** (`src/accelsearch.c`, `src/accel_utils.c`, `src/accel_gather.c`, `src/responses.c`, `src/corr_routines.c`, `src/maximize_rz.c`, `src/maximize_rzw.c`, `include/accel.h`)
- Fourier-domain acceleration search (Ransom, Eikenberry & Middleditch 2002), with an optional jerk search (Andersen & Ransom 2018) and harmonic summing of 1, 2, 4, 8, 16 or 32 harmonics. — [docs/accelsearch.1](B/docs/accelsearch.1)
- Options and defaults: `-numharm 8` (1–32), `-zmax 200` (0–1200), `-wmax` (0–4000, off by default), `-sigma 2.0`, `-flo 1.0` Hz, `-fhi 10000` Hz (bounds apply to the highest harmonic), `-rlo/-rhi` (honoured since v6.3.0), `-inmem` (whole f-fdot plane in RAM), `-photon`, `-median`, `-locpow`, `-nonorm` (v6.3.3), `-zaplist`, `-baryv`, `-otheropt`, `-noharmpolish`, `-noharmremove`, `-ncpus`, `-perfinfo`. — [clig/accelsearch_cmd.cli](B/clig/accelsearch_cmd.cli)
- Outputs: `<root>_ACCEL_<zmax>` (text), `<root>_ACCEL_<zmax>.cand` (binary `fourierprops`), `.txtcand`. With `-wmax` the names are `_ACCEL_<z>_JERK_<w>`. The `.inf` is appended to the text file (`cat %s.inf >> accelnm`). — [src/accel_utils.c](B/src/accel_utils.c)
- Text columns: Cand, Sigma, Summed Power, Coherent Power, Num Harm, Period (ms), Frequency (Hz), FFT 'r', Freq Deriv (Hz/s), FFT 'z', FFT 'w', Accel (m/s^2), Notes. A summary of fundamentals comes first, then per-harmonic details. — [src/accel_utils.c](B/src/accel_utils.c) `output_fundamentals`
- `.cand` record = `fourierprops` {double r; float rerr; double z; float zerr; double w; float werr; float pow, powerr, sig, rawpow, phs, phserr, cen, cenerr, pur, purerr, locpow}. — [include/presto.h](B/include/presto.h)
- Reported f, fdot (and fdotdot) are averages over the observation, i.e. at the midpoint. prepfold's reference is the start. — [FAQ.md](B/FAQ.md)
- `.dat` input is allowed: accelsearch FFTs and dereddens it itself, and `-zaplist`/`-baryv` zap it. — [docs/accelsearch.1](B/docs/accelsearch.1)
- Tutorial guidance: 8 harmonics is almost always enough for MSPs, zmax should stay under about 200, and run time doubles per harmonic stage and grows linearly with zmax. — [tutorial](B/docs/PRESTO_search_tutorial.pdf). For a 5–15% orbit-fraction regime use jerk. — [Andersen & Ransom 2018, arXiv:1807.07900](https://arxiv.org/abs/1807.07900)
- GPU forks (not part of PRESTO): PrestoZL (recommended in the FAQ, includes jerk), presto_on_gpu, Laidler's fork, and AstroAccelerate FDAS. — [FAQ.md](B/FAQ.md)

**ACCEL_sift.py / sifting.py** (`examplescripts/ACCEL_sift.py`, `python/presto/sifting.py`)
- Globs `*ACCEL_*0` and `*DM*.inf`, then applies `read_candidates` → `remove_duplicate_candidates` → `remove_DM_problems(min_num_DMs=2, low_DM_cutoff=2.0)` → `remove_harmonics` → writes the list sorted by sigma. Script settings: `sigma_threshold 4.0` (module default 6.0), `c_pow_threshold 100`, `r_err 1.1` bins, `short_period 0.0005` s, `long_period 15` s, `harm_pow_cutoff 8.0`, and `known_birds_p/_f`. — [examplescripts/ACCEL_sift.py](B/examplescripts/ACCEL_sift.py), [python/presto/sifting.py](B/python/presto/sifting.py)
- Harmonic removal checks integer ratios 1–16 in both directions plus 3/2, 5/2, 2/3, 4/3, 5/3, 3/4, 5/4, 2/5, 3/5 and 4/5, with `f_err = r_err/T`. — [sifting.py](B/python/presto/sifting.py)
- The DM comes from the filename via `DM_re = r"_DM(\d+\.\d{2,3})(?!\d)"` (2–3 decimals since v6.1.0). Output is `candfile:candnum` for prepfold, with ASCII sigma-vs-DM plots. — [sifting.py](B/python/presto/sifting.py), [tutorial](B/docs/PRESTO_search_tutorial.pdf). ACCEL_sift sigmas are single-trial so that they can be compared across DMs. — [FAQ.md](B/FAQ.md)

**prepfold / prepfold_multi / show_pfd** (`src/prepfold.c`, `src/prepfold_pipeline.c`, `src/prepfold_utils.c`, `src/prepfold_plot.c`, `src/fold.c`, `src/polycos.c`, `src/prepfold_multi.c`, `src/show_pfd.c`, `python/presto/prepfold.py`)
- Folds raw data, `.dat`, `.subNN` files or events. Candidates are specified by `-p/-pd/-pdd`, `-f/-fd/-fdd`, `-accelcand N -accelfile X.cand` (`-pfact`/`-ffact` for harmonics), `-par`, `-psr`, `-timing <par>` (tempo2 polycos, no search), `-polycos`, or `-bin` orbit. It writes `.pfd`, `.pfd.bestprof`, `.pfd.ps` and `.png` (Ghostscript). — [docs/prepfold.1](B/docs/prepfold.1), [clig/prepfold_cmd.cli](B/clig/prepfold_cmd.cli)
- Search defaults: `-npart 64`, `-pstep 2`, `-pdstep 4`, `-dmstep 2` (profile bins), `-npfact 2` (2·npfact·proflen+1 period and pdot trials), `-ndmfact 3` (2·ndmfact·proflen+1 DMs). Flags: `-nosearch`, `-nopsearch`, `-nopdsearch`, `-nodmsearch`, `-searchpdd`, `-searchfdd`, `-fine`, `-coarse`, `-slow`, `-fixchi`, `-mask`, `-ignorechan`. `-n` defaults to the samples per period. — [clig/prepfold_cmd.cli](B/clig/prepfold_cmd.cli)
- Trial spacing is set by phase drift `Δφ = Δf T + Δḟ T²/2`, quantized in profile bins, with about ±2 wraps for p/pd and about ±3 rotations across the band for DM. DM never goes negative. — [FAQ.md](B/FAQ.md)
- Since v6.3.0 raw-data folds longer than 2 h (or with `-usepoly`) use tempo2 polycos instead of the cubic `bary2topo()` fit. The phase error of the cubic fit is about 0.002 turns at 2 h, 0.1 at 8 h and 3 turns at 20 h for a 2 ms pulsar. — [CHANGELOG.md](B/CHANGELOG.md) v6.3.0
- `prepfold_multi` (v6.0.0) folds many candidates in one raw-data pass with OpenMP, and with polycos since v6.3.0. — [CHANGELOG.md](B/CHANGELOG.md)
- `.pfd` binary header: ints `numdms, numperiods, numpdots, nsub, npart, proflen, numchan, pstep, pdstep, dmstep, ndmfact, npfact`; length-prefixed strings `filenm, candnm, telescope, pgdev`; 16-char `rastr`/`decstr`; doubles `dt, startT, endT, tepoch, bepoch, avgvoverc, lofreq, chan_wid, ...`. The trial arrays and the profile cube plus `foldstats` follow. — [src/prepfold_utils.c](B/src/prepfold_utils.c) `write_prepfoldinfo`. `show_pfd` re-plots and zaps a `.pfd` (RFI removal after the fact). — [tutorial](B/docs/PRESTO_search_tutorial.pdf), [FAQ.md](B/FAQ.md). Python `presto.prepfold.pfd` reads it. — [FAQ.md](B/FAQ.md)

**single_pulse_search.py** (`bin/single_pulse_search.py`)
- Boxcar matched filtering of `.dat` files through FFT convolution, with no phase dependence. Boxcars are `[1, 2, 3, 4, 6, 9, 14, 20, 30, 45, 70, 100, 150, 220, 300]` bins. The default uses ≤30 (`max_downfact = 30`); `-m/--maxwidth` (s) enables the larger ones. Kernels are normalized by `1/sqrt(downfact)`. `chunklen = 8000`. Detrending uses blocks of `-d` × 1000 samples (piecewise linear by `scipy.signal.detrend`; `-f/--fast` uses median removal for about 2× speed). The block std is a trimmed std (central 95%) × 1.148. Bad blocks have stds outside median ± 4σ of the std distribution (`-b` disables this), and since v6.1.2 dead blocks (std ≤ 1e-3 median) are skipped. The default threshold is `-t 5.0`. Duplicates are pruned within downfact/2. — [bin/single_pulse_search.py](B/bin/single_pulse_search.py)
- Definition: `sigma = sum(signal - bkgd)/RMS/sqrt(boxcar_width)`, which stays the same under downsampling while the pulse is resolved. — [bin/single_pulse_search.py](B/bin/single_pulse_search.py) docstring
- `.singlepulse` ASCII columns: `# DM  Sigma  Time (s)  Sample  Downfact`. Re-running on `.singlepulse` files only re-plots. — [bin/single_pulse_search.py](B/bin/single_pulse_search.py)
- v6.1.0 searches the tail samples that used to be dropped (up to 7999 per file). — [CHANGELOG.md](B/CHANGELOG.md). The `presto.singlepulse` subpackage (make_spd.py, plot_spd.py, rrattrap.py) was removed in v6.1.0. — [CHANGELOG.md](B/CHANGELOG.md)

**explorefft / exploredat** (`src/explorefft.c`, `src/exploredat.c`)
- Interactive PGPLOT viewers of `.fft` and `.dat` files, used to find birdie widths. A `.zaplist` can be loaded into explorefft. — [tutorial](B/docs/PRESTO_search_tutorial.pdf). explorefft had an under-reported power at high |z| until v6.1.0. — [CHANGELOG.md](B/CHANGELOG.md)

**get_TOAs.py** (`bin/get_TOAs.py`, `python/presto/fftfit.py`, `fftfit_aarchiba.py`)
- Works on `.pfd` files folded with `-timing`/par or `-nosearch` (or `-u` since v6.3.0 for searched folds). Template: `-g FWHM` gaussian, a `.bestprof`, or a `pygaussfit.py` multi-gaussian. `-n` TOAs must factor npart, and `-s` gives subband TOAs. — [tutorial](B/docs/PRESTO_search_tutorial.pdf), [CHANGELOG.md](B/CHANGELOG.md). FFTFIT is pure Python since v6.0 (`-A classic|aarchiba`), and CHIME support was added in v6.0.1. — [CHANGELOG.md](B/CHANGELOG.md)

**coherentsearch / coherent_sift.py** (v6.2.0; `src/coherentsearch.c`, `src/cs_interp.c`, `src/cs_chunk.c`, `src/cs_metric.c`, `src/cs_search.c`, `src/cs_plot.c`, `bin/coherent_sift.py`)
- Coherent harmonic summing over normalized `.fft` files (run `rednoise` first). It is a C port of CoherentSearch.jl. It has no acceleration search. — [docs/coherentsearch.1](B/docs/coherentsearch.1), [CHANGELOG.md](B/CHANGELOG.md)
- Defaults: `-nharms 60`, `-m 16` (interpolation kernel bins), `-hidr 0.5` (trial spacing at the highest harmonic), `-threshold 8.0` S/N, `-maxdecim 6`, `-lofreq 0.1`, `-hifreq 125` Hz, `-ncands 100`, `-drtol 1.0`, `-numharm 16` (for harmonic collapse), `-blocksize 2048`, boxcar width spacing factor 1.5 up to 0.3 of the profile. `-deep` sets `-nharms 120 -maxdecim 12 -hifreq 5` → `<input>_deep.cohout`. — [clig/coherentsearch_cmd.cli](B/clig/coherentsearch_cmd.cli), [src/cs_interp.c](B/src/cs_interp.c) `cs_default_params`
- Output is `<input>.cohout` (CoherentSearch.jl format) with S/N, freq, period, harmonics and duty cycle. Known pulsars are labelled. `-plot` gives `.ps`/`.png`. Comb-like candidates (strongest harmonic ≥ 7, with >8× the power of harmonics 1–4) are dropped unless `-noperiodic`. — [docs/coherentsearch.1](B/docs/coherentsearch.1)

**search_bin** (`src/search_bin.c`, `src/minifft.c`, `src/bincand.c`)
- Phase-modulation ("sideband") binary search (Ransom, Cordes & Eikenberry 2003, ApJ 589, 911). It runs many power-of-two "mini-FFTs" (`-minfft 32` to `-maxfft 65536`, `-overlap 0.25`) of harmonically summed (`-harmsum 3`) powers, with `-flo 1` to `-fhi 2000` Hz and `-ncand 100`. Candidates are refined with `bincand`. The man page states that it "has not (yet!) led to the discovery of a new pulsar". — [docs/search_bin.1](B/docs/search_bin.1), [clig/search_bin_cmd.cli](B/clig/search_bin_cmd.cli)

**FFA: absent**
- No FFA program or source exists in `src/`, `bin/` or `python/presto/`. A grep for "ffa" hits only unrelated identifiers such as `-ffact`. ROADMAP.md mentions FFA only as riptide, an external comparison for coherentsearch. — [ROADMAP.md](B/ROADMAP.md); repo tree at commit 4a5bd80

**Other utilities seen in `bin/`**
- `stacksearch.py` (stacked FFT search, v5.1), `fourier_fold.py` (profiles from `.fft` complex amplitudes, v5.1), `pfdzap.py`, `quick_prune_cands.py`, `rfifind_stats.py` + `weights_to_ignorechan.py` (build `-ignorechan` lists), `waterfaller.py`, `injectpsr.py`, `psrfits2fil.py`, `downsample_filterbank.py`, `detrend_dat.py`, `pygaussfit.py`, `sum_profiles.py`, `fitorb.py`, `compare_periods.py`. — repo [bin/](https://github.com/scottransom/presto/tree/master/bin), [CHANGELOG.md](B/CHANGELOG.md)

### Inferences
- File naming carries metadata, e.g. the DM in `_DM12.34` and zmax in `_ACCEL_200`. Downstream sifting parses filenames (regex), so renaming files breaks the pipeline.
- `accelsearch` defaults (sigma 2, numharm 8, zmax 200) are deliberately permissive. Real thresholding happens in sifting (sigma 4–6) and folding.

### Gaps
- I did not read `mpiprepsubband.c` in detail (the MPI broadcast scheme comes from the tutorial and FAQ only).
- I did not read `explorefft.c`/`exploredat.c` key bindings or `show_pfd` options.
- I did not read `zapbirds` interactive mode details.

---

## 2. Standard workflows (flowchart-ready) with file extensions

### Takeaway
The tutorial gives a 13-step canonical search: readfile → rfifind → prepdata DM0 → realfft/explorefft → accelsearch -zmax 0 for birdies → .birds/zaplist → DDplan.py → prepsubband → realfft (+ rednoise) + zap → accelsearch (or coherentsearch) → ACCEL_sift.py → prepfold → get_TOAs.py. single_pulse_search.py branches off the `.dat` files.

### Cited Findings
- Tutorial outline: "1) Examine data format (readfile) 2) Search for RFI (rfifind) ... 4) FFT the time series (realfft) 5) Identify 'birdies' ... 6) Make zaplist (makezaplist.py ...) 7) Make De-dispersion plan (DDplan.py) 8) De-disperse (prepsubband) 9) Search ... (accelsearch) 10) ... single pulses (single_pulse_search.py) 11) Sift ... (ACCEL_sift.py) 12) Fold ... (prepfold) 13) Start timing ... (prepfold and get_TOAs.py)". — [PRESTO_search_tutorial.pdf](B/docs/PRESTO_search_tutorial.pdf)
- Example commands from the tutorial: `rfifind -time 2.0 -o Lband GBT_Lband_PSR.fil`; `prepdata -nobary -o Lband_topo_DM0.00 -dm 0.0 -mask Lband_rfifind.mask ...`; `realfft Lband_topo_DM0.00.dat`; `accelsearch -numharm 4 -zmax 0 Lband_topo_DM0.00.dat`; `makezaplist.py Lband.birds`; `zapbirds -zap -zapfile Lband.zaplist -baryv -5.697334e-05 *.fft` (or `simple_zapbirds.py Lband.birds *.fft`); `prepsubband -nsub 32 -lodm 0.0 -dmstep 2.0 -numdms 24 -downsamp 4 -mask Lband_rfifind.mask -o Lband GBT_Lband_PSR.fil`; `accelsearch -zmax 0 Lband_DM0.00.fft`; `python ACCEL_sift.py > cands.txt`; `prepfold -accelcand 1 -accelfile Lband_DM62.00_ACCEL_0.cand Lband_DM62.00.dat`; `prepfold -n 64 -nsub 96 -p 0.004621638 -dm 62.0 GBT_Lband_PSR.fil`; `single_pulse_search.py *.dat`; `get_TOAs.py -g 0.1 -n 20 newpulsar.pfd`; coherent: `rednoise` → `coherentsearch -plot *_red.fft` → `coherent_sift.py --html cands.html . > cohcands.txt`. — [tutorial](B/docs/PRESTO_search_tutorial.pdf)
- Periodicity pipeline (file flow):
  raw `.fil`/`.fits` → **rfifind** → `_rfifind.mask` (+ `.stats`, `.bytemask`, `.rfi`, `.inf`, `.ps`)
  raw + mask → **prepsubband** (plan from **DDplan.py**) → `<root>_DM<dm>.dat` + `.inf` (optional `.subNN`)
  `.dat` → **realfft** → `.fft` → (**rednoise** → `_red.fft` + `_red.inf`) → (**zapbirds**/**simple_zapbirds.py** with `.birds` → `.zaplist`, in place)
  `.fft` → **accelsearch** → `_ACCEL_<z>` (text) + `_ACCEL_<z>.cand` (binary) [+ `_JERK_<w>`]
  all ACCEL → **ACCEL_sift.py** → ranked `cands.txt` (`candfile:candnum`)
  cand + `.dat` / `.subNN` / raw → **prepfold** → `.pfd`, `.pfd.bestprof`, `.pfd.ps`, `.png`
  Alternative branch: `_red.fft` → **coherentsearch** → `.cohout` (+ `_deep.cohout`) → **coherent_sift.py** → prepfold `-f <freq> -dm <DM>`. — assembled from [tutorial](B/docs/PRESTO_search_tutorial.pdf), [src/accel_utils.c](B/src/accel_utils.c), [src/rfifind.c](B/src/rfifind.c), [src/prepsubband.c](B/src/prepsubband.c), [docs/coherentsearch.1](B/docs/coherentsearch.1)
- Single-pulse pipeline: `.dat` (+ `.inf`) → **single_pulse_search.py** → `<root>_DM<dm>.singlepulse` + `_singlepulse.ps`. Re-plot with `single_pulse_search.py *DM1??.??*.singlepulse`. `waterfaller.py` makes raw-data waterfalls. DDplan downsampling should keep dispersion smearing ≲ 1 bin. — [tutorial](B/docs/PRESTO_search_tutorial.pdf), [bin/single_pulse_search.py](B/bin/single_pulse_search.py)
- Folding/timing pipeline: candidate → `prepfold` (search) → confirm on raw data with `-nsub 64+` → make a par → `prepfold -timing X.par raw` (polycos via tempo2, no search) → `.pfd` → `get_TOAs.py -g FWHM -n N` → TOAs (TEMPO/tempo2 format). Folding `.dat` is fast for triage, and folding `.subNN` near the candidate's DM is faster than raw. — [tutorial](B/docs/PRESTO_search_tutorial.pdf), [docs/prepfold.1](B/docs/prepfold.1)
- Example scripts: `examplescripts/dedisp.py`, `full_analysis.py`, `ACCEL_sift.py`, `short_analysis_simple.py`, `ppdot_plane_plot.py`; survey scripts `bin/GBNCC_search.py`, `PALFA_presto_search.py`, `GBT350_drift_search.py`. — [CHANGELOG.md](B/CHANGELOG.md) v3.0; [bin/](https://github.com/scottransom/presto/tree/master/bin)
- Hint: running `-ncpus` does not speed up most tools. Parallelize by running one process per DM, and set `OPENBLAS_NUM_THREADS=1` when running many jobs. — [FAQ.md](B/FAQ.md); [CHANGELOG.md](B/CHANGELOG.md) v6.2.0

### Inferences
- The tutorial's two-stage approach is to make birdies on the DM 0 topocentric series and then zap with `-baryv` on the barycentred DM trials. It exists because birdies are topocentric while searched series are barycentric.

### Gaps
- The `tests/test_realdata_search.py` pipeline was not read (it mirrors the tutorial according to the CHANGELOG).

---

## 3. Algorithms and math inside the programs

### Takeaway
PRESTO relies on frequency-domain statistics throughout. Powers are normalized to unit mean noise (|A|²/median·ln2 or by local blocks). The f-ḟ(-f̈) plane is built by correlating short FFT segments with precomputed Fresnel-integral template kernels at Δr = 0.5 (interbinning), Δz = 2, Δw = 20. Harmonic summing adds subharmonic planes at rounded fractional (r, z, w). Significance is converted from a χ² with 2·numharm dof, with a trials correction, to an "equivalent Gaussian sigma". prepfold uses epoch-folding reduced χ² with an empirically corrected effective dof.

### Cited Findings
**Power normalization / rednoise**
- Block-median normalization in the search: `norm = 1/sqrt(median(|A|^2)/ln 2)` applied to amplitudes over each fftlen/2 data block. For exponential (χ²₂) noise the median equals ln 2 × the mean. Other options are photon `1/sqrt(nph)` (freq-0 normalization), `-locpow` double-tophat local means, and `-nonorm`. — [src/accel_utils.c](B/src/accel_utils.c) `subharm_fderivs_vol`
- v6.3.3: `_red.fft` (and `.dat`) input that passes a noise-scale check (within 10%) is searched as `-nonorm`. Programs stop with an error if the scale is off by more than 2× (`check_fft_normalization`). — [CHANGELOG.md](B/CHANGELOG.md) v6.3.3
- rednoise `dered_engine`: bin 0 is set to 1+0i. Block lengths grow logarithmically, `bufflen = startwidth*ln(binnum)` capped at `endwidth`, until `binnum/T ≥ endfreq` Hz. Each block's level is `median(power)/ln 2`, and amplitudes are divided by `sqrt(level)` linearly interpolated (`dslope`) between neighbouring block midpoints. — [src/misc_utils.c](B/src/misc_utils.c)

**Fourier interpolation / interbinning / responses**
- Interbinning: `ACCEL_NUMBETWEEN 2`, `ACCEL_DR 0.5` (half-bin steps). The data are "spread" (`spread_no_pad`) and correlated so that half-integer bins are computed. — [include/accel.h](B/include/accel.h), [src/accel_utils.c](B/src/accel_utils.c)
- Kernel half-widths: `z_resp_halfwidth = int(0.5*1.1*|z|) + NUMFINTBINS` (LOWACC) or `+3*NUMFINTBINS` (HIGHACC), with `NUMFINTBINS 16`. The w (jerk) half-width uses the parabola extremes `r0 = 0.5(w/6 − z)`, `r1 = 0.5(w/6 + z)`. Responses are built from Fresnel integrals (`fresnl`). — [src/responses.c](B/src/responses.c), [include/presto.h](B/include/presto.h)
- Definitions: r = fT (Fourier bin), z = ḟT² (bins drifted over the observation), w = f̈T³. — [include/presto.h](B/include/presto.h) `fourierprops` comments; paper: [Ransom, Eikenberry & Middleditch 2002, arXiv:astro-ph/0204349](https://arxiv.org/abs/astro-ph/0204349) (AJ 124, 1788: Fourier interpolation, recovering power lost to frequency drift, estimating ḟ, f̈, centroid, purity, significance)
- Correlation search: `ACCEL_USELEN 7470` good points per 8K FFT (works to zmax 300). Since v6.3.3 the FFT lengths come from a measured `next_good_fftlen()` list (powers of 2 plus 7-smooth lengths), and `-nonorm` with small kernels uses 1280-point fundamental FFTs. — [include/accel.h](B/include/accel.h), [CHANGELOG.md](B/CHANGELOG.md)
- Grid: `ACCEL_DZ 2`, `ACCEL_DW 20`. Independent candidates must be farther apart than `ACCEL_CLOSEST_R 15.0` bins. — [include/accel.h](B/include/accel.h)

**Harmonic summing**
- Stages are numharm = 2, 4, 8, 16, 32. At stage N the planes of harmonics (2k+1)/N, k = 0..N/2−1, are added to the running fundamental sum (the even fractions were added in earlier stages). The required sub-plane coordinates are `r_sub = rint(2·r·h/N)·0.5`, `z_sub = nearest(0.5·z·h/N)·2`, `w_sub = nearest(0.05·w·h/N)·20`. Each stage is thresholded at `powcut[stage]`. — [src/accel_utils.c](B/src/accel_utils.c) `calc_required_r/z/w`, `add_search_ffdotpows`
- `-zmax` refers to the highest harmonic. Subharmonics drift proportionally less. `-flo/-fhi` refer to the highest harmonic, so candidates can fall outside them. — [tutorial](B/docs/PRESTO_search_tutorial.pdf), [FAQ.md](B/FAQ.md)
- Candidate optimization (`optimize_accelcand`): per-harmonic or joint ("harmpolish", default) maximization `max_rz_arr_harmonics` / `max_rzw_arr_harmonics` (GSL Nelder–Mead since v6.0, replacing NR amoeba). The bottom-of-file powers use local power that excludes the 5 closest bins on each side. — [src/accel_utils.c](B/src/accel_utils.c), [src/maximize_rz.c](B/src/maximize_rz.c), [FAQ.md](B/FAQ.md), [CHANGELOG.md](B/CHANGELOG.md)

**Sigma, trials, thresholds** (`src/characteristics.c`)
- `candidate_sigma(power, numsum, numtrials)`: χ² = 2·power with dof = 2·numsum. `logp = chi2_logp(χ², dof) + ln(numtrials)` is converted by `equivalent_gaussian_sigma(logp) = gsl_cdf_gaussian_Qinv(exp(logp))`, with an asymptotic extension for logp < −600 and a floor of 0. — [src/characteristics.c](B/src/characteristics.c)
- `chi2_logp` uses `gsl_cdf_chisq_Q` unless χ² > 1200 or (dof > 30 and χ² > 100), in which case it uses the asymptotic incomplete-gamma expansion (A&S 26.4.19). The transitions were improved in v6.0.0. — [src/characteristics.c](B/src/characteristics.c), [CHANGELOG.md](B/CHANGELOG.md)
- `power_for_sigma(sigma, numsum, numtrials) = 0.5·gsl_cdf_chisq_Qinv(Q_gauss(sigma)/numtrials, 2·numsum)`. — [src/characteristics.c](B/src/characteristics.c)
- Number of independent trials per stage (numharm = 2^i): zmax = 0 gives `(rhi−rlo)/numharm`. Acceleration gives `(rhi−rlo)·(numz+1)·(dz/6.95)/numharm`. Jerk gives `(rhi−rlo)·(numz+1)·(dz/6.95)·(numw+1)·(dw/44.2)/numharm`. — [src/accel_utils.c](B/src/accel_utils.c) `create_accelobs`
- FAQ: top-of-file sigmas include this single-run trials correction, while bottom (per-harmonic) sigmas are single-trial. The default `-sigma 2` "you probably don't want to go much below 1.5". The rough spectral S/N per harmonic is sqrt(P) − 1 (not recommended). — [FAQ.md](B/FAQ.md)
- rfifind power threshold: `pow_reject = power_for_sigma(freqsig, 1, ptsperint/2)`. Mean and std outliers are flagged when |value − interval median| or |value − channel median| > `timesig × global robust std`. The robust statistics come from `calc_avgmedstd` with a 0.8 (avg, std) or 0.5 (pow) kept fraction, and interval/channel medians that deviate by > timesig·σ are replaced by the global median. Then whole-interval masking applies when bad channels > chanfrac·numchan, and whole-channel masking when bad intervals > intfrac·numint. — [src/rfifind_plot.c](B/src/rfifind_plot.c)

**prepfold folding and χ²**
- `fold()` treats each sample as finite-duration and spreads it over the profile bins it covers. It supports f, ḟ, f̈, delays (binary) and on/off pairs. — [src/fold.c](B/src/fold.c)
- Per-profile statistics: `prof_var = data_var·numdata·profbinwidth` and `redchi = Σ(prof − prof_avg)²/(prof_var·(numprof − 1))`. — [src/fold.c](B/src/fold.c)
- Effective dof: `DOF_corr(dt_per_bin) = dt_per_bin·0.96·(1 + dt_per_bin^1.806)^(−1/1.806)` (Monte Carlo fit after Demorest/Brisken), where dt_per_bin = P/nbins/dt. The reported significance is `chi2_logp(redchi·dofeff, dofeff)` → equivalent sigma, labelled "Prob(Noise) < ... (~X sigma)" in `.bestprof`. — [src/fold.c](B/src/fold.c), [src/prepfold_plot.c](B/src/prepfold_plot.c)
- The search does not refold. It shifts the npart × nsub sub-profiles by phase delays (period: `ii·pstep·Δ/npart`; pdot: `fdot2phasedelay`; fdotdot; DM: subband delays via `dm_from_delay(dphase·dmstep, f)`) and combines them (`combine_profs`), maximizing reduced χ². `-fixchi` rescales so the off-pulse reduced χ² is 1 (`estimate_offpulse_redchi2`). — [src/prepfold_pipeline.c](B/src/prepfold_pipeline.c), [src/prepfold_utils.c](B/src/prepfold_utils.c)
- FAQ: the χ² significance is single-trial and valid only if the off-pulse reduced χ² is about 1 (epoch folding, Leahy et al. 1983). — [FAQ.md](B/FAQ.md)

**Dedispersion**
- `delay = DM/(0.000241·f²)` s. Channel frequencies are channel centres, Doppler-corrected by the average v/c when barycentring. Subband two-stage dedispersion is described in Section 1. — [src/dispersion.c](B/src/dispersion.c), [src/prepdata.c](B/src/prepdata.c)

**Coherent harmonic summing (coherentsearch)**
- Harmonic h of trial t sits at `h·r_lo + h·t·pnum/q` bins. It is interpolated exactly with an m = 16-bin kernel of weights `1/(dr − j)` and scale `conj(sin(π dr) e^{iπ dr})/π` (RE&M 2002 Eq. 30), with the weights tabulated by residue. The harmonic stack (60 by default) is inverse-FFT'd to a profile. Score: `max_{w,p} (Σ_{i=p}^{p+w−1} P_i − (w/nbins) Σ P_i)·shape_w/σ` over a geometric boxcar bank (×1.5 steps up to 0.3 of the profile, wrapping). σ is computed analytically from which harmonics are filled, assuming unit-mean powers (or measured by MAD with `-measuredsigma`). The decimation ladder folds every k-th harmonic at k× the fundamental. — [src/cs_interp.c](B/src/cs_interp.c), [src/cs_metric.c](B/src/cs_metric.c), [src/cs_chunk.c](B/src/cs_chunk.c), [ROADMAP.md](B/ROADMAP.md)
- Paper: "A Coherent Harmonic Summing Pulsar Search Code", S. M. Ransom, arXiv:2610.11749, submitted 2026-10-08. In over 860,000 injections at matched false-alarm rate in white noise it recovers 76%, against 64% for riptide FFA (recommended configuration) and 42% for accelsearch. Its threshold is stable under red noise. It is 1.4–2.2× faster than riptide on one core, about 27× faster on 48 cores, and up to 9× faster on GPUs than a 20-core workstation (GPU path in the Julia code). — [arXiv:2610.11749](https://arxiv.org/abs/2610.11749). ROADMAP.md quotes riptide at "50–71%" and accelsearch `-zmax 0` at 42%. — [ROADMAP.md](B/ROADMAP.md) (a different summary of the riptide configurations; not a contradiction of the 64% "recommended" figure, but note the range)

**Jerk search**
- Andersen & Ransom 2018 (ApJL 863, L13; arXiv:1807.07900) extends FDAS to constant jerk. It improves sensitivity when T ≈ 5–15% of P_orb at much higher cost, and it found Ter5am. — [arXiv:1807.07900](https://arxiv.org/abs/1807.07900). Added in PRESTO v2.1. — [CHANGELOG.md](B/CHANGELOG.md)

### Inferences
- The factors 6.95 and 44.2 in the trials count act as effective correlation lengths (in z and w units) of the response, i.e. independent cells are about 6.95 z units and 44.2 w units wide. The source does not explain these values.

### Gaps
- I did not read the exact derivation of the χ²₂-to-σ extension (`extended_equiv_gaussian_sigma`) beyond its call.
- I did not read the Ransom 2001 PhD thesis. The 2002 AJ paper text beyond its abstract was not fetched, so claims about template correlation rely on the code and the man page.

---

## 4. Version changes 2024–2026 and license

### Takeaway
PRESTO moved to meson in v5.0 (Apr 2024) and became an ERFA/tempo2-based, Fortran-free, conda-forge package in v6.0 (Jul 2026). It added coherentsearch in v6.2 (Nov 2026), polyco folding in v6.3.0, and three accelsearch/coherentsearch speed releases in v6.3.1–6.3.3 (14–17 Nov 2026). The license is GPL-2.0.

### Cited Findings
- Release dates (GitHub releases API): v5.0.0 2024-04-21; v5.0.1 2024-05-21; v5.0.2 2024-06-22; v5.0.3 2025-02-13; v5.1.0 2025-08-27; v5.2.0 2025-10-29; v5.3.0/5.3.1 2026-03-05; v6.0.0 2026-07-05; v6.0.1 2026-08-18; v6.1.0 2026-09-10; v6.1.1 2026-10-14; v6.1.2 2026-11-09; v6.2.0 2026-11-10; v6.3.0 2026-11-14; v6.3.1 2026-11-14; v6.3.2 2026-11-15; v6.3.3 2026-11-17. — [GitHub releases](https://github.com/scottransom/presto/releases)
- v5.0.0: meson + meson-python build (numpy.distutils deprecated), Python ≥ 3.8, Spigot code and slalib removed. — [CHANGELOG.md](B/CHANGELOG.md)
- v5.0.2: Numpy 2.0 C API and Python ≥ 3.9. v5.1.0: stacksearch.py, fourier_fold.py, pfdzap.py, DDplan smearing fix. v5.2.0: GSL required, and K&R declarations updated for GCC 15 (with AI-assisted tooling, prompted by Paul Ray). v5.3.0: rfifind OpenMP, accelsearch speed-ups, and fixes for >2³² / >2³¹ point data. — [CHANGELOG.md](B/CHANGELOG.md)
- v6.0.0: no `$PRESTO` environment variable needed (data in `{prefix}/share/presto`); TEMPO no longer needed (ERFA barycentring; polycos via `tempo2 -tempo1 -polyco`); Fortran-free (FFTFIT in Python, GSL least squares, GSL Brent); SLALIB and RANDLIB replaced by ERFA and GSL; median via GSL (even-length median changed); NR amoeba replaced by GSL Nelder–Mead; new `prepfold_multi`; pyclig man pages; fast `makewisdom` by default (`-patient` for the old behaviour); `dftfold` removed; conda-forge package `presto-pulsar`. — [CHANGELOG.md](B/CHANGELOG.md)
- v6.1.0: zapped bins filled with noise by default; f-fdot optimizer kernel-width bug (dating from 2001) fixed, with recovered power at |z| = 200 rising from 0.5995 to 0.9994, although the search sigmas were unaffected; deredden moved to libpresto; `rednoise -inplace`, `realfft -rednoise`; SIGPROC fix for declinations between 0 and −1°; single_pulse_search tail fix; presto.singlepulse removed. — [CHANGELOG.md](B/CHANGELOG.md)
- v6.1.1: `prepfold -accelcand` jerk with `-start` fixed (`f += fd·dt + fdd·dt²/2`); signed 8-bit on aarch64 fixed. v6.1.2: accelsearch about 2.4× faster on non-FMA builds (out-of-line `fmaf()` removed); `mask.c` out-of-bounds read that could change which channels were masked; single_pulse dead blocks; 1/2/4-bit `presto.filterbank`; CI with sanitizers and weekly real-data tests (Terzan 5). — [CHANGELOG.md](B/CHANGELOG.md)
- v6.2.0: coherentsearch and coherent_sift.py; Python ≥ 3.12; pkgconf; long plot filenames. — [CHANGELOG.md](B/CHANGELOG.md)
- v6.3.0: polyco folding (`-usepoly`, automatic > 2 h), prepfold_multi with par files, `get_TOAs.py -u`, accelsearch on FFTs with ≥ 2³¹ bins (PR #233), `-rlo/-rhi` honoured (#232), known-pulsar matching fixes. v6.3.1: runtime SIMD dispatch (AVX2/AVX-512) for coherentsearch, 1.2–1.35× faster; `-perfinfo`. — [CHANGELOG.md](B/CHANGELOG.md)
- v6.3.2: accelsearch threads over spectrum blocks (8 threads about 5.5× faster than 1; output byte-identical for any `-ncpus`); vectorized threshold scan (1.2–1.45×); `-ftree-vectorize`. — [CHANGELOG.md](B/CHANGELOG.md)
- v6.3.3: with one thread about 1.25× faster than 6.3.2, about 1.6× faster than v5.0.0, and about 8× faster than v5.3.1–v6.1.1 (8 threads 7.7× faster than v5.0.0), from fused harmonic-sum stages (`add_search_ffdotpows`) and an AVX2/AVX-512/NEON gather. Also measured FFT lengths (rerun `makewisdom`; `makewisdom -check`), `-nonorm`, and automatic no-block-median for rednoised input (6–11% fewer false candidates above σ = 2 on noise). Benchmark: `-zmax 200 -numharm 8` on a 100M-bin spectrum, 16-core Sapphire Rapids. — [CHANGELOG.md](B/CHANGELOG.md), [docs/accel_speed_releases.png](B/docs/accel_speed_releases.png)
- Unreleased after v6.3.3: prepfold fix for short `.dat` files, RA/Dec "60 seconds" rounding fix, ppgplot exceptions, SIGPROC missing-position memory corruption (#234). — [CHANGELOG.md](B/CHANGELOG.md)
- License: GNU GPL version 2 (`COPYING`); CITATION.cff says `license: GPL-2.0`; the GitHub API reports `GPL-2.0`. Since v6.1.0 the conda package is "plain GPL-2.0-only" (ERFA is no longer bundled). ASCL entry: ascl.net/1107.017. — [COPYING](B/COPYING), [CITATION.cff](B/CITATION.cff), [CHANGELOG.md](B/CHANGELOG.md)
- Install: `conda install -c conda-forge presto-pulsar` (linux-64, osx-64, linux-aarch64, osx-arm64). — [CHANGELOG.md](B/CHANGELOG.md) v6.0–6.1.2

### Inferences
- The burst of releases in Nov 2026 (6.1.2 through 6.3.3 within 9 days) indicates very active development. Results from different versions are not byte-comparable (zapping noise from v6.1.0, rednoised normalization from v6.3.3, even-length median from v6.0).

### Gaps
- There are no release notes for v5.0.3's exact date context beyond the API. Release dates are GitHub publication times, not tag times.

---

## 5. Why PRESTO's intermediate files are hard to reuse elsewhere

### Takeaway
The formats are headerless or custom binary and depend on sidecar files, filename conventions and in-code normalization conventions. They follow PRESTO's in-memory C structs and carry no versioning or self-description, so other tools must re-implement PRESTO-specific conventions to read them correctly.

### Cited Findings
- `.dat` is raw float32 with no header. All metadata lives in a separate ASCII `.inf` (the `infodata` struct: N, dt, MJD int+frac, RA/Dec, DM, low channel freq, total BW, channel width, num_chan, bary flag, on/off pairs, telescope, notes). — [include/makeinf.h](B/include/makeinf.h); FAQ entry on reading `.dat`/`.fft` into Python — [FAQ.md](B/FAQ.md)
- `.fft` stores only N/2 complex64 values, with Nyquist packed into the imaginary part of bin 0. Its length and dt come from `.inf`, and `_red.fft` changes the normalization (DC set to 1). — [src/fastffts.c](B/src/fastffts.c), [src/misc_utils.c](B/src/misc_utils.c)
- `.mask`, `.stats`, `.bytemask`, `.pfd` and `.cand` are raw dumps of C struct fields in native byte order. For example, `.pfd` readers have byte-swap code (`swap_double` in `read_prepfoldinfo`), and `.cand` is an array of `fourierprops`. — [src/mask.c](B/src/mask.c), [src/rfifind.c](B/src/rfifind.c), [src/prepfold_utils.c](B/src/prepfold_utils.c), [include/presto.h](B/include/presto.h)
- Semantics are embedded in the code. Examples: z and w are in Fourier units (ḟT², f̈T³); accelsearch f/ḟ refer to the midpoint, while prepfold refers to the start; the `.inf` epoch is the top-channel centre at the first sample; sigma means trial-corrected at the top of the file and single-trial at the bottom. — [FAQ.md](B/FAQ.md)
- Pipelines parse metadata from filenames: the `_DM\d+\.\d{2,3}` regex, `_ACCEL_<z>` suffixes, and `.subNN` numbering. — [python/presto/sifting.py](B/python/presto/sifting.py), [src/prepsubband.c](B/src/prepsubband.c)
- `.dat` files are by default barycentred with samples added/dropped (not resampled), so they are not uniformly topocentric. Birdie zapping needs the observation's average v/c (`-baryv`). — [src/prepdata.c](B/src/prepdata.c), [tutorial](B/docs/PRESTO_search_tutorial.pdf)
- The block size for SIGPROC data is hard-coded at 2400 spectra, which constrains downsampling factors and rfifind interval granularity. — [FAQ.md](B/FAQ.md)
- coherentsearch's analytic noise "assumes exactly" PRESTO's `realfft` → `rednoise` normalization. — [ROADMAP.md](B/ROADMAP.md)

### Inferences
- Interoperating with other ecosystems (e.g. SIGPROC `.tim`, PSRCHIVE archives, YOUR/sigpyproc) requires converters: PRESTO ships `dat2tim.py`, `tim2dat.py`, `psrfits2fil.py` and `pfd_for_timing.py` (all in `bin/`). Its Python readers (`presto.infodata`, `presto.prepfold`, `presto.rfifind`) are the practical reference implementations.
- Because results changed across 6.x (zapping, normalization), stored intermediates should record the PRESTO version. The files themselves do not.

### Gaps
- I did not verify whether `.inf` contains a version field (the struct shown has none). I did not inspect `dat2tim.py`/`tim2dat.py` details.
