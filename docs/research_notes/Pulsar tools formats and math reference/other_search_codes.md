# Non-PRESTO pulsar / FRB search codes: internals, pipelines, formats and math (as of Nov 2026)

Method note: almost everything below was read directly from shallow clones of each code's original repository (HEAD as of 2026-11-24), plus the papers' arXiv abstract pages. Line-level claims cite the repo file. "Last commit" dates are from `git log -1` on the default branch or from the GitHub API `pushed_at` field.

## Q1. Exact step order, inputs/outputs, header fields, core algorithms and math of each code

### Takeaway
Every single-pulse code here follows the same skeleton: read and clean the data, dedisperse, remove the baseline and normalise, run a boxcar matched filter, threshold, then cluster or merge duplicates. They differ in the dedispersion engine (brute-force GPU dedisp in heimdall and peasoup, AA's GPU direct kernels, CPU subband in TransientX, FDMT in dmt and PulsarX), in the duplicate-merging method (heimdall's friends-of-friends labelling vs TransientX's DBSCAN), and in whether an ML stage comes after (FETCH, ASTROFLOW and DRAFTS). On the periodicity side, peasoup is an FFT code (time-domain resampling for acceleration, then FFT and harmonic sum), riptide is an FFA code, and PulsarX is the folder and optimiser that sits downstream of both.

### Cited Findings

#### 1. heimdall (Barsdell, Jameson et al.; GPU single-pulse search)
- **Repository and status.** The original code is on SourceForge (`git.code.sf.net/p/heimdall-astro/code`). Its last commit is dated 2023-06-15 ("fix incorrect power calculation in sigproc 32fp unpacker"), so the code is in low-activity maintenance. No canonical `ajameson/heimdall-astro` GitHub repo exists (the API returns 404). — [heimdall SourceForge repo](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/)
- **License.** It is ambiguous. `COPYING` contains the GNU GPL text, but source headers (e.g. `Pipeline/hd/are_coincident.cuh`) say "Copyright (C) 2012 by Ben Barsdell and Andrew Jameson. Licensed under the Academic Free License version 2.1". — [are_coincident.cuh](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Pipeline/hd/are_coincident.cuh)
- **Hardware and dependencies.** It needs a CUDA GPU, Thrust, and the `dedisp` library (Barsdell's GPU brute-force direct dedispersion). PSRDADA ring-buffer input is optional. — [heimdall configure/config/*.m4](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/config/)
- **Inputs.** SIGPROC filterbank (`-f`) or a PSRDADA ring buffer (`-k`). From the SIGPROC header, `Formats/SigprocFile.C` reads `nchans`, `nbits`, `ibeam`, `tsamp` (s), `fch1`, `foff` and `tstart` (MJD, converted to a UTC `time_t`). — [Formats/SigprocFile.C](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Formats/SigprocFile.C)
- **CLI options** (`Pipeline/parse_command_line.C`): `-dm lo hi`, `-dm_tol`, `-dm_pulse_width`, `-dm_nbits`, `-no_scrunching`, `-scrunch_tol`, `-rfi_tol`, `-rfi_min_beams`, `-rfi_no_narrow`, `-rfi_no_broad`, `-boxcar_max`, `-detect_thresh`, `-cand_sep_time`, `-cand_sep_filter`, `-cand_sep_dm_trial`, `-cand_rfi_dm_cut`, `-max_giant_rate`, `-min_tscrunch_width`, `-baseline_length`, `-nsamps_gulp`, `-zap_chans`, `-boxcar_renorm`, `-coincidencer host:port`, `-fswap`, `-gpu_id`, `-beam`, `-output_dir`. — [parse_command_line.C](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Pipeline/parse_command_line.C)
- **Defaults** (`Pipeline/default_params.C`): `nsamps_gulp=262144`, `baseline_length=2.0` s, `dm_tol=1.25`, `dm_pulse_width=40` µs, `dm_nbits=32`, `use_scrunching=true`, `scrunch_tol=1.15`, `rfi_tol=5.0`, `boxcar_max=4096` samples, `detect_thresh=6.0`, `cand_sep_time=3`, `cand_sep_filter=3` (filter index), `cand_sep_dm=200` (DM trials, not DM units), `min_tscrunch_width=4096`, `max_giant_rate=0` (no limit). — [default_params.C](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Pipeline/default_params.C)
- **Pipeline order in `hd_execute()`** (`Pipeline/pipeline.cu`, for each gulp):
  1. `clean_filterbank_rfi()` does zero-DM narrow- and broad-band RFI excision. Masked samples are replaced with values sampled at random from nearby data, which produces the killmask.
  2. `apply_manual_killmasks()` applies the `-zap_chans` list.
  3. `dedisp_set_killmask` and `dedisp_execute_adv` dedisperse every DM trial. The DM list comes from `dedisp_generate_dm_list(dm_min, dm_max, pulse_width, tol)`, with adaptive time scrunching (`dedisp_enable_adaptive_dt`), so each DM has a scrunch factor. `nsamps_computed = nsamps − max_delay`.
  4. A per-DM loop copies the series to the GPU.
  5. `RemoveBaselinePlan.exec` removes the baseline, with smoothing radius `baseline_length/(2·dt)`.
  6. `GetRMSPlan` normalises the series to unit RMS.
  7. `MatchedFilterPlan.prep` takes a prefix sum, then a loop over boxcar widths `filter_width = scrunch, 2·scrunch, … ≤ boxcar_max` (powers of 2) applies the boxcar. Above `min_tscrunch_width` it decimates in time. Each filtered output is normalised by 1/√w, or renormalised by measured RMS if `-boxcar_renorm`.
  8. `GiantFinder.exec(thresh=detect_thresh)` finds giants. If `max_giant_rate` is exceeded, the search aborts.
  9. After all DMs, `label_candidate_clusters()` and then `merge_candidates()` group the giants.
  10. Output is written as text or sent to the coincidencer socket.

  — [pipeline.cu](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Pipeline/pipeline.cu)
- **Baseline removal math** (`remove_baseline.cu`). The series is resampled to a grid whose length is `sample_count·5^k`. It is then repeatedly median-of-5 scrunched down to `sample_count ≈ 2·N/(2·smooth_radius)` points (oversample=2), and the ends are linearly extrapolated. The result is linearly stretched back to N samples as the baseline and subtracted. This gives a robust running-median-like detrend, and the code comments say it is "robust against outliers (e.g., strong RFI spikes)". — [remove_baseline.cu](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Pipeline/remove_baseline.cu)
- **Matched filter math** (`matched_filter.cu`). An inclusive prefix scan is taken once, S[k] = Σ_{j<k} x[j]. The boxcar of width w is then y[i] = S[i+ahead] − S[i−behind], with ahead=⌈w/2⌉ and behind=⌊w/2⌋. Time scrunching is done by striding through S. So S/N = y/(σ√w), since σ=1 after normalisation. — [matched_filter.cu](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Pipeline/matched_filter.cu)
- **Giant finding** (`find_giants.cu`). It extracts all samples above `thresh` and splits them into isolated "giants" using a `merge_dist` (head flags from `not_nearby`, then an inclusive scan for segment ids). It then computes each giant's peak S/N, index, begin and end with `reduce_by_key`. — [find_giants.cu](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Pipeline/find_giants.cu)
- **Candidate clustering math** (`hd/are_coincident.cuh`). Two giants i and j are coincident if |samp_i−samp_j| ≤ time_tol·max(2^filter_i, 2^filter_j), |filter_i−filter_j| ≤ filter_tol and |dm_i−dm_j| ≤ dm_tol (in trial indices). `label_candidate_clusters.cu` resolves equivalence chains on the GPU (each label is replaced by the minimum label of the root chain) to give friends-of-friends connected components. `merge_candidates.cu` sorts by label and uses `reduce_by_key`, keeping the max-S/N member and counting members. — [are_coincident.cuh](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Pipeline/hd/are_coincident.cuh), [label_candidate_clusters.cu](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Pipeline/label_candidate_clusters.cu), [merge_candidates.cu](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Pipeline/merge_candidates.cu)
- **Output format.** One text file per gulp, `<output_dir>/<UTC YYYY-MM-DD-hh:mm:ss>_<beam two-digit>.cand`. It has 9 tab-separated columns: `snr`, `sample_idx`, `time_s` (= sample_idx·dt), `filter_idx` (boxcar width = 2^filter_idx samples), `dm_trial_idx`, `dm`, `members`, `begin_sample`, `end_sample`. The same rows go to the coincidencer socket, prefixed by a header line. Coincidencer tools are in `Applications/coincidencer.C` and `Applications/candidate_profiler.C`, and plotting scripts in `Scripts/trans_*.py`. — [pipeline.cu](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Pipeline/pipeline.cu)
- **Paper.** Barsdell et al. 2012, MNRAS 422, 379, "Accelerating incoherent dedispersion" (the dedisp library). Heimdall itself has no dedicated paper. — [arXiv:1201.5380](https://arxiv.org/abs/1201.5380)
- **Python reimplementation.** sigpyproc 2.0 ships `algorithms/heimdall.py`, a "standalone, inspectable reference implementation" of heimdall's boxcar, giant-finding and clustering stages. — [sigpyproc heimdall.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/algorithms/heimdall.py)

#### 2. AstroAccelerate (Armour group, Oxford; GPU library and standalone)
- **Status and license.** GPL-3.0. GitHub `pushed_at` is 2026-09-25 and the last default-branch commit is 2025-12-01. The code is CUDA-only (nvcc) and builds with Make, CMake or Meson. It ships a C++ API (`aa_pipeline_api`) and a Python wrapper (`python/`). — [AstroAccelerate repo](https://github.com/AstroAccelerateOrg/astro-accelerate), [MANUAL.md](https://github.com/AstroAccelerateOrg/astro-accelerate/blob/master/MANUAL.md)
- **Components.** The components are `dedispersion` (output: a dedispersed time chunk as `std::vector<unsigned short>`), `analysis` (single-pulse search), `periodicity` (FFT, harmonic sum and peak find) and `fdas` (Fourier-domain acceleration search), plus a jerk search. The API validates the component combination and fixes the execution order. The user binds a *plan* per component, and AA converts it into a GPU-optimised *strategy* (e.g. `aa_ddtr_strategy`). — [MANUAL.md](https://github.com/AstroAccelerateOrg/astro-accelerate/blob/master/MANUAL.md)
- **Configuration file** (`input_files/ska_test_file.txt`):
  - `range lo hi step inBin outBin` lines give the piecewise DDTR plan with time binning.
  - Other keys: `sigma_cutoff`, `sigma_constant`, `max_boxcar_width_in_sec`, `periodicity_sigma_cutoff`, `periodicity_harmonics`, `selected_card_id`.
  - Flags, enabled by removing the leading "-": `analysis`, `acceleration`, `acceleration_jerk`, `z_max/z_step/w_max/w_step`, `periodicity`, `output_dmt`, `zero_dm`, `zero_dm_with_outliers`, `rfi`, `oldrfi`, `threshold`, `baselinenoise`, `fdas_custom_fft`, `fdas_inbin`, `fdas_norm`, `set_bandpass_average`, `copy_ddtr_data_to_host`, and `file <path.fil>`.

  — [ska_test_file.txt](https://github.com/AstroAccelerateOrg/astro-accelerate/blob/master/input_files/ska_test_file.txt)
- **Execution order** (`include/aa_permitted_pipelines_generic.hpp`, `run_pipeline()`): `set_data` → `input_data_preprocessing` (zero-DM, zero-DM-with-outliers and RFI per flags; bandpass normalisation) → DDTR per time chunk and DM range (`aa_device_dedispersion_kernel.cu`, with `aa_corner_turn`, `aa_bin_gpu`) → SPD analysis per chunk → after all chunks, `periodicity()`, `acceleration()` (FDAS) and `jerk_search()`. — [aa_permitted_pipelines_generic.hpp](https://github.com/AstroAccelerateOrg/astro-accelerate/blob/master/include/aa_permitted_pipelines_generic.hpp)
- **Header fields read** (`src/aa_sigproc_input.cpp`): `tstart`, `tsamp`, `nbits`, `fch1`, `foff`, `nchans`, `nifs` (the reader rejects `nifs != 1`), and others. nsamp is computed from the file size and nbits (32, 16 or 8 bits and lower). A PRESTO-input reader (`aa_presto_input.cpp`) also exists. — [aa_sigproc_input.cpp](https://github.com/AstroAccelerateOrg/astro-accelerate/blob/master/src/aa_sigproc_input.cpp)
- **SPDT (single-pulse detection) math.**
  - `SPDT_search_long_MSD_plane()` in `src/aa_device_SPS_long.cu` builds a boxcar-filter bank iteratively. The first plane computes boxcars up to `nBoxcars` taps.
  - Each later iteration decimates the series by 2 (`d_decimated`) and reuses the partial boxcar sums (`d_boxcar_values`), starting from `startTaps`. Widths therefore grow while each output sample is evaluated at decimated resolution.
  - The S/N is computed with mean and σ interpolated from an MSD (mean/standard-deviation) "plane profile" (`d_MSD_interpolated`) that depends on boxcar width. These statistics come from `aa_device_MSD*.cu`, which can optionally use outlier rejection (`aa_device_MSD_outlier_rejection_kernel.cu`; the `baselinenoise` flag).
  - The maximum S/N over widths and its tap count are stored per (DM, t), then thresholding and peak-finding follow (`aa_device_threshold.cu`, `aa_device_peak_find.cu`).
  - The algorithm and its sensitivity trade-offs (boxcar vs IIR, etc.) are published in Adámek & Armour 2020, ApJS 247, 56.

  — [aa_device_SPS_long.cu](https://github.com/AstroAccelerateOrg/astro-accelerate/blob/master/src/aa_device_SPS_long.cu), [README](https://github.com/AstroAccelerateOrg/astro-accelerate/blob/master/README.md)
- **SPS output.** Binary float32 quadruplets `(DM, time_s = TS·tsamp + tstart, SNR, width = BW·inBin)` written to `analysed-t_<tstart>-dm_<lo>-<hi>.dat`. With peak-finding enabled the file is `peak_analysed-t_…dat`. — [aa_device_analysis.cu](https://github.com/AstroAccelerateOrg/astro-accelerate/blob/master/src/aa_device_analysis.cu)
- **FDAS.** This is the Fourier-domain correlation technique (as in Ransom's accelsearch) on the GPU: the spectrum is convolved with a bank of z (and optionally w, jerk) templates by overlap-save cuFFT, which gives the f–ḟ plane. Code is in `aa_fdas_device.cu`, `aa_fdas_host.cu`, `aa_device_acceleration_fdas.cu` and `aa_device_jerk_search.cu`. Outputs are text `acc_list_<DM>.dat` (candidate list) and `acc_<DM>.dat` (full f-fdot plane if `output_ffdot_plan`). The reference paper is Dimoudi et al. 2018, ApJS 239, 28. — [aa_fdas_host.cu](https://github.com/AstroAccelerateOrg/astro-accelerate/blob/master/src/aa_fdas_host.cu), [README](https://github.com/AstroAccelerateOrg/astro-accelerate/blob/master/README.md)
- **Dedispersion paper.** Novotný et al. 2023, ApJS 269, 1 ("Accelerating Dedispersion Using Many-core Architectures"). — [README](https://github.com/AstroAccelerateOrg/astro-accelerate/blob/master/README.md)

#### 3. dmt (pravirkr/dmt, a Dispersion Measure Transform library)
- **Status and license.** MIT. Very active: the last commit is 2026-11-20, v0.8.0 was released 2026-11-11 and v0.7.0 on 2026-11-03. It is C++20 with Python bindings (`dmtlib`) and needs Python ≥3.12 and CUDA ≥12.2. — [dmt README](https://github.com/pravirkr/dmt), [CHANGELOG](https://github.com/pravirkr/dmt/blob/main/CHANGELOG.md)
- **Engines.** These share one plan, DM-grid and convention layer:
  - FDMT (tree) and FDMT-FFT
  - DDMT (direct brute force), SDMT, DDMT-FFT
  - CohFDMT (coherent)

  The CPU path uses OpenMP and SIMD (AVX2, AVX-512, NEON). The GPU path uses CUDA and HIP and gives "the same results to float rounding". Streaming blocks are joined exactly through an internal history. The code reads packed 1–16-bit input directly. Every plan gives "the exact sigma of each DM trial, so S/N needs no calibration run". — [dmt README](https://github.com/pravirkr/dmt)
- **Source files.** `lib/fdmt.cpp`, `lib/fdmt_plan*.cpp`, `lib/fdmt_fft*.cpp`, `lib/ddmt*.cpp`, `lib/sdmt.cpp`, `lib/cfdmt*.cpp`, `lib/boxcar_*.cpp`, `lib/fdmt_boxcar.cpp`, `lib/kadane_plan.cpp` and `lib/dm_grids.cpp`. Python bindings are in `src/dmt_python.cpp`, and benchmarks in `bench/*_b.cpp/.cu`. — [dmt repo tree](https://github.com/pravirkr/dmt)
- **API conventions.** The waterfall is `(nchans, nsamps)` in **ascending** frequency. For SIGPROC files with foff<0, use `dmtlib.sigproc_band(fch1, foff, nchans)` and `dmtlib.to_ascending()`. By default, output times refer to `f_min`. `compute_fdmt(...)` returns the DM–time plane plus a plan, and `dmtlib.dedisperse(..., engine="fdmt_fft")` gives white-noise time series for periodicity searches. The C++ class is `dmt::algorithms::FDMT(f_min, f_max, nchans, nsamps, tsamp, dt_max, dt_min, dt_step, use_box_smearing, mode, Exec::cpu|cuda|hip)`. The dispersion constant is the TEMPO2 value 1/2.41e-4. — [dmt README](https://github.com/pravirkr/dmt), [CHANGELOG](https://github.com/pravirkr/dmt/blob/main/CHANGELOG.md)
- **Benchmark preset.** The `meerkat_lband` preset is 4096 channels, 856–1712 MHz, 612.48 µs sampling, 20.1 s blocks, and up to 4057 DM trials over 0–2000 pc cm⁻³, run on an EPYC 9654 and an L40S. — [dmt README](https://github.com/pravirkr/dmt)
- **FDMT math** (Zackay & Ofek 2017, ApJ 835, 11). The algorithm works up a tree. It starts from single channels with delay range Δt ∈ [0, Δt_max(sub-band)]. At each level, adjacent sub-band pairs [f0,f1] and [f1,f2] are merged:

  A_{[f0,f2]}(Δt, t) = A_{[f0,f1]}(Δt₁, t) + A_{[f1,f2]}(Δt − Δt₁, t − Δt₁),

  where Δt₁ is the part of the delay that falls inside the first half for that DM curve, Δt ∝ DM(f⁻² − f'⁻²). There are log₂N_f iterations, and the complexity is about N_t·N_Δ·log₂N_f (versus N_t·N_f·N_Δ for brute force). — [Zackay & Ofek 2017, arXiv:1411.5373](https://arxiv.org/abs/1411.5373)

#### 4. TransientX (Men & Barr; CPU single-pulse search)
- **Status and license.** GPL-3.0. The last commit is 2026-07-23 and the code is actively maintained. It is C++ (autotools) and needs Boost, PlotX, XLibs and SOFA/ERFA. A Docker image is available (`ypmen/pulsarx`). — [TransientX repo](https://github.com/ypmen/TransientX)
- **Paper.** Men & Barr 2024, A&A 683, A183 ("TransientX: A high performance single pulse search package"). It describes a CPU-only real-time search with sub-band dedispersion, kd-tree DBSCAN clustering, about 9% S/N loss for a boxcar against a scattered Gaussian, and a note on spectral-index mismatch for multi-octave receivers. — [arXiv:2401.13834](https://arxiv.org/abs/2401.13834)
- **Inputs.** SIGPROC filterbank or PSRFITS (`--psrfits`, plus `--wts`, `--scloffs` and `--zero_off` for PSRFITS scaling), read in chunks of `-l` seconds. — [docs/pipeline.rst](https://github.com/ypmen/TransientX/blob/master/docs/pipeline.rst)
- **Pipeline order** (`SinglePulse::prepare/run` in `src/singlepulse/singlepulse.cpp`):
  1. `downsample` (`--td`, `--fd`)
  2. `equalize`: per-channel zero mean and unit variance
  3. `baseline`: subtract a smoothed zero-DM profile (`--baseline` width)
  4. `rfi`: a list of filters run in order, `zap flow fhigh` (manual), skewness–kurtosis IQR filter, `mask td fd` (IQR outlier mask), `kadaneF td fd`, `kadaneT`, `zero` and `zdot`
  5. `dedisp`: subband dedispersion driven by a ddplan
  6. `boxcar`: matched filter
  7. `cluster`: DBSCAN
  8. `candplot`: PNG plus candidate list. `dedisp.rundump` can optionally write dedispersed series.

  — [singlepulse.cpp](https://github.com/ypmen/TransientX/blob/master/src/singlepulse/singlepulse.cpp), [docs/pipeline/*.rst](https://github.com/ypmen/TransientX/tree/master/docs/pipeline)
- **ddplan format.** Columns are `td fd dms ddm ndm snrloss maxwidth(s)`, for example `1 1 0.0 0.1 3344 0.1 0.1` and `2 1 334.4 0.2 1255 0.1 0.1`. — [docs/pipeline/dedispersion.rst](https://github.com/ypmen/TransientX/blob/master/docs/pipeline/dedispersion.rst)
- **RFI filters.**
  - **zdot** (Men et al. 2019, MNRAS 488, 3957) subtracts, per channel, the component correlated with the zero-DM profile: x_c(t) ← x_c(t) − (⟨x_c·z⟩/⟨z·z⟩)·z(t).
  - The **skewness–kurtosis** filter (Men et al. 2023) removes channels whose skewness or kurtosis are IQR outliers.
  - **kadaneF** finds the maximum-sum window along time per channel (Kadane's algorithm), to remove slow-varying RFI.
  - **mask** is an IQR outlier mask in the time–frequency plane.

  — [docs/pipeline/rfi_mitigation.rst](https://github.com/ypmen/TransientX/blob/master/docs/pipeline/rfi_mitigation.rst)
- **Boxcar math** (`boxcar.cpp`). Widths run from `--minw` to `--maxw` (seconds). Successive widths grow geometrically by `wfactor = 1/(1−snrloss)²`, so neighbouring widths lose at most a fraction `snrloss` of S/N, and each step is at least one sample. S/N uses a `sqrt(1/(wn·var))` normalisation. `--iqr` gives a robust noise estimate. — [boxcar.cpp](https://github.com/ypmen/TransientX/blob/master/src/singlepulse/boxcar.cpp)
- **Clustering** (`cluster.cpp`). Above-threshold (`--thre`) points in (time, DM-smearing) space go into a 2-D kd-tree and then DBSCAN. `eps` is `-r` (ms, with a minimum of 1.5·tsamp) and `min_samples` is `-k`. `--maxncand` keeps the top-S/N clusters per chunk, and `--drop` removes candidates whose width equals maxw (likely RFI). — [cluster.cpp](https://github.com/ypmen/TransientX/blob/master/src/singlepulse/cluster.cpp), [docs/pipeline/clustering.rst](https://github.com/ypmen/TransientX/blob/master/docs/pipeline/clustering.rst)
- **Output.** A `.cands` text file with columns `ibeam`, `cand_no`, `mjd` (15 dp), `dm`, `width_ms`, `snr`, `fl`, `fh` (frequency extent in MHz), `png_name`, `id`, `filename`. There is also one PNG per candidate, and with `--saveimage`, images for AI classifiers. `replot_fil` re-reads `.cands`, re-folds, applies cuts (`--dmcutoff`, `--widthcutoff`, `--snrcutoff`), removes duplicates and can write PSRFITS archives (`archivewriter_tx.cpp`). — [candplot.cpp](https://github.com/ypmen/TransientX/blob/master/src/singlepulse/candplot.cpp), [README](https://github.com/ypmen/TransientX)
- **Usage example.** `transientx_fil -v -o test -t 4 --zapthre 3.0 --fd 1 --overlap 0.1 --ddplan ddplan.txt --thre 7 --maxw 0.1 --snrloss 0.1 -l 2.0 --drop -z kadaneF 8 4 zdot -f *.fil`. JSON config is also accepted (`-c`). — [README](https://github.com/ypmen/TransientX)

#### 5. PulsarX (Men et al.; CPU folding, candidate optimisation and RFI/dedispersion toolkit)
- **Status and license.** GPL-3.0. The last commit is 2026-11-08 and the code is actively maintained. It is C++ and shares the XLibs and PlotX dependencies with TransientX. — [PulsarX repo](https://github.com/ypmen/PulsarX)
- **Paper.** Men et al. 2023, A&A 679, A20. Its two key ideas are the **pruned FDMT (pFDMT)**, a "cache-friendly" FDMT that "reuses intermediate processing results and prunes the unused computation paths", and a **Tikhonov-regularised least-squares folding method (TLSM)** that improves profile time resolution. For about 500 candidates, dedispersion operations drop by a factor of about 50 compared with brute force. The code is used in the MeerKAT MMGPS and TRAPUM surveys. — [arXiv:2309.02544](https://arxiv.org/abs/2309.02544)
- **Programs** (`src/pulsar/*.cpp`):
  - `filtool`: cleans a filterbank or PSRFITS file and writes a new filterbank, driven by a `filplan.json`.
  - `dedisperse_all_fil`: RFI cleaning plus a ddplan, then writes dedispersed series (`--format presto` gives `.dat/.inf`).
  - `psrfold_fil`, `psrfold_fil2`, `psrfold_tim`: fold and optimise candidates.
  - `dmffdot`: re-optimise and plot archives in DM/F0/F1.
  - `pulsarsearch`, `gridsearch`: the DM/F0/F1/F2 grid optimiser.
  - Also `pulsar_injector`, `filmaker`, `psrheader`, `archivewriter`, `archivelite`, `dedispersionlite(U)`.

  — [PulsarX src/pulsar](https://github.com/ypmen/PulsarX/tree/main/src/pulsar)
- **psrfold_fil inputs.**
  - Single candidate: `--dm`, `--f0`, `--f1`, `--f2`, `--acc`, `--pepoch`.
  - Binary orbit: `--binary` with `--PB`, `--A1`, `--T0`, `--ECC`, `--OM`.
  - Batch input: `--parfile` or `--t2pred`, or `--candfile` with columns `#id dm acc F0 F1 S/N`.
  - Folding layout: `--template`, `-b/--nbin`, `--nbinplan`, `-L/--tsubint`, `-n/--nsubband`, `--blocksize`.
  - Cleaning: CLFD (`--clfd`, `--clfd2`, `--bandcorr`), RFI (`-z mask|kadaneF|kadaneT|zdot`), patching (`--threPatch`, `--widthPatch`, `--fillPatch`).
  - Search control: `--nosearch`, `--nodmsearch`, `--nof0search`, `--nof1search`, `--f2search`, `--scale`.
  - Outputs: PSRFITS fold archives (`.ar`), PNG plots and an updated `<rootname>_<date>_<beam>.cands` file.

  — [psrfold_fil.cpp](https://github.com/ypmen/PulsarX/blob/main/src/pulsar/psrfold_fil.cpp), [examples/example.cands](https://github.com/ypmen/PulsarX/blob/main/examples/example.cands)
- **Typical use.** `psrfold_fil` is the folder that peasoup's README recommends for peasoup candidates. — [peasoup README](https://github.com/ewanbarr/peasoup)

#### 6. riptide (Morello; FFA periodicity search)
- **Status and license.** MIT. The last commit is 2026-09-20, the latest tag is v0.2.7, and the code is maintained. Its core is C++ with pybind11 and its pipeline is in Python. It runs on the CPU only, parallelised with one process per DM trial. — [riptide repo](https://github.com/v-morello/riptide)
- **Paper.** Morello et al. 2020, MNRAS 497, 4654, "Optimal periodicity searching: revisiting the Fast Folding Algorithm for large-scale pulsar surveys". It shows that the FFA combined with an optimal boxcar S/N reaches near-optimal sensitivity, especially for long periods and small duty cycles where FFT plus incoherent harmonic summing loses sensitivity. — [arXiv:2004.03701](https://arxiv.org/abs/2004.03701)
- **Inputs.** Dedispersed time series in PRESTO format (`.inf/.dat`) or SIGPROC format (`.tim`). The SIGPROC reader parses keys such as `src_raj`, `tstart`, `tsamp`, `nbits`, `fch1` and `refdm` (`reading/sigproc.py`). For SIGPROC `.tim` files, `fmin`, `fmax` and `nchans` must be given in the YAML config because the header lacks band information. They are used to build an internal DM plan that selects a minimal subset of DM trials. — [config/example.yaml](https://github.com/v-morello/riptide/blob/master/src/riptide/pipeline/config/example.yaml), [reading/sigproc.py](https://github.com/v-morello/riptide/blob/master/src/riptide/reading/sigproc.py)
- **Pipeline order** (`rffa` app; `pipeline/pipeline.py` `Pipeline.process`):
  1. `prepare`: DM selection, applying `dmselect.min/max` and `dmsinb_max` (a cap on DM·|sin b|, default 45).
  2. `search`: per DM trial, dereddening with a running-median subtraction (`rmed_width` s, `rmed_minpts`), normalisation, then `ffa_search` per period `range`, then peak detection on the periodogram.
  3. `cluster_peaks`: clustering by frequency.
  4. `flag_harmonics`.
  5. `apply_candidate_filters`.
  6. `build_candidates`: re-fold, and per-candidate the best DM.
  7. `save_products`: `peaks.csv`, `clusters.csv`, `candidates.csv`, `candidate_NNNN.json` and optional PNG plots.

  — [pipeline.py](https://github.com/v-morello/riptide/blob/master/src/riptide/pipeline/pipeline.py)
- **`ffa_search` parameters.** `period_min=1.0`, `period_max=30.0`, `fpmin=8` (minimum number of signal periods), `bins_min=240`, `bins_max=260`, `ducy_max=0.20`, `wtsp=1.5` (geometric step of trial boxcar widths), `deredden=True`, `rmed_width=4.0`, `rmed_minpts=101`. — [search.py](https://github.com/v-morello/riptide/blob/master/src/riptide/search.py)
- **Periodogram and "ranging" math** (`cpp/periodogram.hpp`). The FFA folds at integer bins-per-period. To cover a continuous period range, the series is first downsampled by `ds_ini` such that `ds_ini·tsamp·bins_min = period_min`. For each downsampling factor, the FFA is run for every integer b in [bins_min, bins_max]. The downsampling factor then grows geometrically, with the number of cycles chosen to reach `period_max`. `ceilshift(rows, cols, pmax) = ceil(cols·(rows−1)·(1−cols/pmax))` limits which FFA rows are evaluated. — [periodogram.hpp](https://github.com/v-morello/riptide/blob/master/src/riptide/cpp/periodogram.hpp)
- **FFA transform** (`cpp/transforms.hpp`). The series is reshaped into m rows of p bins. Pairs of halves are merged recursively, with each row of the second half circularly shifted by the appropriate integer, so that all m trial periods in [p, p+1) samples are folded in O(m·p·log₂m) operations rather than O(m²·p). — [transforms.hpp](https://github.com/v-morello/riptide/blob/master/src/riptide/cpp/transforms.hpp), [Morello et al. 2020](https://arxiv.org/abs/2004.03701)
- **S/N math** (`cpp/snr.hpp`, `snr1`). For a folded profile of n bins and trial width w, the boxcar template has zero mean and unit square sum, with height h = √((n−w)/(n·w)) and baseline −b, where b = (w/(n−w))·h. Using a circular prefix sum, dmax = max over phase of the w-bin window sum, and S/N = ((h+b)·dmax − b·Σx)/σ_noise. This is an exact matched-filter S/N for a top-hat pulse on a flat baseline. — [snr.hpp](https://github.com/v-morello/riptide/blob/master/src/riptide/cpp/snr.hpp)
- **Python API.** `TimeSeries` (`normalise`, `deredden`, `downsample`, `fold`, `generate`, `from_binary`), `ffa_search`, `ffa2`, `ffafreq`, `Periodogram` and `Candidate`. — [time_series.py](https://github.com/v-morello/riptide/blob/master/src/riptide/time_series.py)

#### 7. peasoup (Barr et al.; GPU FFT acceleration search)
- **Status and license.** Apache-2.0. The last commit is 2025-08-05, so the code is maintained at a low cadence. It needs CUDA (up to 12.6 according to its README), cuFFT and dedisp. It does **not** fold: the README says "use … `psrfold_fil` (recommended), `prepfold`, or `dspsr`". — [peasoup README](https://github.com/ewanbarr/peasoup)
- **Newer features** (README): template-bank coherent Keplerian searches (circular orbits with three parameters; elliptical orbits with the BT or ELL8 resamplers); segmented searches (`--start_sample`, `--fft_size`, `--nsamples`); `--cdm` (finer acceleration steps near a coherent DM); `--single_precision_harmonic_sums`; dedispersion-only output of PRESTO-compatible `.dat` files (`-d --nosearch`); and a progress bar (`-p`). — [peasoup README](https://github.com/ewanbarr/peasoup)
- **Pipeline order** (`src/pipeline_multi.cu`):
  1. Read the SIGPROC filterbank (`include/data_types/filterbank.hpp`, `header.hpp`).
  2. Build the DM list with `Dedisperser.generate_dm_list(dm_start, dm_end, dm_pulse_width, dm_tol)` and apply an optional killmask.
  3. GPU dedispersion over gulps (`Dedisperser.dedisperse`).
  4. A `DMDispenser` hands DM trials to per-GPU `Worker`s. Per DM, the worker does an R2C cuFFT, forms the power spectrum, takes a running median, dereddens (`Dereddener`), zaps birdies (`Zapper`, from a zaplist), forms an interpolated spectrum, computes stats (mean and rms), then a C2R iFFT back to a whitened time series.
  5. For each acceleration in `AccelerationPlan.generate_accel_list(dm, cdm)`: time-domain resample (`TimeDomainResampler.resampleII`), R2C FFT, interbinned power spectrum (`form_interpolated`), `stats::normalise` (using the mean and std from step 4), harmonic summing (`HarmonicFolder.fold`), `PeakFinder.find_candidates` on the fundamental and each harmonic sum, then `HarmonicDistiller.distill`.
  6. Per DM, `AccelerationDistiller` (or `Template_Bank_Circular_Distiller`).
  7. Globally, `DMDistiller` then `HarmonicDistiller`.
  8. Optional `MultiFolder`/`FoldOptimiser` (`folding_kernels.cu`).
  9. `OutputFileWriter` writes `<outdir>/overview.xml` with search parameters and candidates. The binary `candidates.peasoup` write is commented out in this file.

  — [pipeline_multi.cu](https://github.com/ewanbarr/peasoup/blob/master/src/pipeline_multi.cu)
- **Acceleration-resampling math** (`src/kernels.cu`). `accel_fact = a·tsamp/(2c)`, and the output sample is `out[i] = in[round(i + i·accel_fact·(i − N))]` (`getAcceleratedIndexII`). This removes the quadratic Doppler phase about the series midpoint, i.e. t' ≈ t + a·t(t−T)/(2c), so that constant-acceleration signals become periodic before the FFT. — [kernels.cu](https://github.com/ewanbarr/peasoup/blob/master/src/kernels.cu)
- **Acceleration step math** (`include/utils/utils.hpp`, `AccelerationPlan`).
  - Per-channel dispersive smear: t_dm = 4.15e-3·|DM−cDM|·(f_bot⁻² − f_top⁻²), with f in GHz and the band of one channel (cfreq ± ch_bw/2).
  - Broadening: t_b = √(t_dm² + w² + tsamp²), where w is the `pulse_width` parameter.
  - Step: δa = 2·t_b·24·c/T_obs²·√(tol²−1). With `tol≈1.1`, this gives about 10% smearing.
  - The accel list includes 0 when the range spans 0.

  — [utils.hpp](https://github.com/ewanbarr/peasoup/blob/master/include/utils/utils.hpp)
- **Key headers.** `include/transforms/{dedisperser,resampler,ffter,spectrumformer,dereddener,birdiezapper,harmonicfolder,peakfinder,distiller,folder,foldoptimiser,scorer,template_bank_reader}.hpp` and `include/data_types/{filterbank,header,timeseries,fourierseries,candidates,folded}.hpp`. Python tools are `tools/peasoup_tools.py` and `peasoup_as_text.py`. — [peasoup repo tree](https://github.com/ewanbarr/peasoup)

#### 8. sigpyproc3 / sigpyproc (FRBs org; Python filterbank toolbox)
- **Status and license.** MIT. v2.0.0 was released 2026-10-17 and is the first PyPI release under the name `sigpyproc`. The `legacy` branch keeps the old API and its C++ bindings. Version 2 is pure Python plus **numba** kernels (it depends on numpy≥2, numba≥0.61, astropy, rocket-fft, h5py, bottleneck, click/rich). It runs on the CPU only. — [HISTORY.rst](https://github.com/FRBs/sigpyproc3/blob/main/HISTORY.rst), [pyproject.toml](https://github.com/FRBs/sigpyproc3/blob/main/pyproject.toml)
- **Readers** (`readers.py`): `FilReader` (SIGPROC), `PFITSReader` (search-mode PSRFITS, read row by row; Stokes I = P0+P1 for coherence data; channel order follows `DAT_FREQ`), `open_filterbank()` and `PulseExtractor`. Header handling is in `header.py` and `io/sigproc.py`. Since v2.0, 2- and 4-bit samples are packed least-significant-bit first, matching dspsr. — [readers.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/readers.py), [HISTORY.rst](https://github.com/FRBs/sigpyproc3/blob/main/HISTORY.rst)
- **`Filterbank` (abstract base, `base.py`) methods.** `read_block`, `read_dedisp_block`, `read_plan` (gulped iteration), `compute_stats`/`chan_stats`, `collapse`, `bandpass`, `dedisperse` (to a TimeSeries), `read_chan`, `invert_freq`, `apply_channel_mask`, `zap_times`, `downsample`, `extract_samps/chans/bands`, `requantize`, `remove_zerodm`, `subband`, `fold(period, dm, accel, nbins=50, nints=32, nbands=32)` (gives a FoldedData cube), `tree_dedisperse`, `dedisperse_plan`, `clean_rfi(method="mad"|"iqrm", threshold=3, freq_mask, …)` (returns a new file and an `RFIMask`), and `view` (GUI). — [base.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/base.py)
- **`TimeSeries` and `FourierSeries`.**
  - `TimeSeries` (`timeseries.py`): `normalise`, `downsample`, `pad`, `deredden`, `fold`, `rfft`, `apply_boxcar`, `resample(accel)`, `correlate`, `search`, `find_pulses` (defaults to a modular PRESTO-inspired single-pulse pipeline in v2.0), `to_dat/to_tim`, `from_dat/from_tim`.
  - `FourierSeries` (`fourierseries.py`): `ifft`, `form_spec(interpolate)`, `deredden`, `recon_prof`, `multiply`, `to_fft/to_spec`, `from_fft/from_spec`.
  - `PowerSpectrum`: `bin2freq`, `bin2period`, `harmonic_fold`, `find_peaks`.

  — [timeseries.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/timeseries.py), [fourierseries.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/fourierseries.py)
- **IQRM math** (`core/rfi.py`, `iqrm_mask`; the method is from Morello et al. 2022). For each per-channel statistic x_i and each lag k in ±1..radius, the code forms lagged differences d_ik = x_i − x_{i+k}. It then flags outliers in each lag's distribution using Tukey IQR fences with threshold in sigma units, and votes across lags. `double_mad_mask` is the alternative. — [core/rfi.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/core/rfi.py)
- **Algorithms and CLI.** `algorithms/{box,sps,heimdall,period,common}.py` contain the boxcar, PRESTO-style and heimdall-style single-pulse pipelines and spectral peak picking. The CLI has `spp_header`, `spp_clean`, `spp_decimate`, `spp_dedisperse`, `spp_fold`, `spp_search`, `spp_bandpass`, `spp_io` (merge/splice/extract), `spp_digifil`, `spp_fake` (FRB populations) and `spp_view`. — [sigpyproc src tree](https://github.com/FRBs/sigpyproc3/tree/main/src/sigpyproc)

#### 9. your (thepetabyteproject; "Your Unified Reader")
- **Status and license.** GPL-3.0. Version 0.6.8; GitHub `pushed_at` is 2026-09-19 and the last default-branch commit 2026-09-19 (default branch `main`). It is Python with optional GPU support (cupy/numba in `utils/gpu.py`). Its paper is Aggarwal et al. 2020 in JOSS. — [your repo](https://github.com/thepetabyteproject/your)
- **`Your` class** (`your/your.py`). `class Your(PsrfitsFile, SigprocFile)` presents PSRFITS and SIGPROC filterbank data through one interface. Properties include `chan_freqs`, `native_tsamp`, `native_foff`, `native_nchans`, `native_nspectra` and `tend`. Methods are `bandpass()`, `get_data(nstart, nsamp, time_decimation_factor, frequency_decimation_factor, pol, npoln)` and `dispersion_delay()`. A unified `Header` class exposes `tsamp`, `nchans`, `foff` and `nspectra`. — [your.py](https://github.com/thepetabyteproject/your/blob/main/your/your.py)
- **`Writer` class** (`your/writer.py`). It converts formats or sub-selects data (channel range, time range, decimation). Its methods are `to_fil`, `to_fits` and `to_dada` (with `setup_dada`/`dada_header`, which writes to a PSRDADA buffer). Format back-ends are `your/formats/{filwriter,fitswriter,psrfits,pysigproc,dada}.py`. — [writer.py](https://github.com/thepetabyteproject/your/blob/main/your/writer.py)
- **Candidate maker** (`Candidate(Your)` in `your/candidate.py`, `bin/your_candmaker.py`).
  - Input is a CSV with columns `file`, `snr`, `width` (log2 of the width; the code uses `2**row["width"]`), `dm`, `label`, `stime`, `chan_mask_path`, `num_files`. This matches heimdall's `filter_idx` convention.
  - Processing order: `get_chunk()` (time window covering the dispersion sweep) → `dedisperse()` (CPU or GPU) → `dmtime(dmsteps=256)` (DM–time butterfly from 0 to 2·DM) → optionally `optimize_dm()` (maximise S/N over DM) → `decimate(key="ft"/"dmt", factor = width//2, mode="median")` → `resize` to 256×256 by default (`--time_size`, `--frequency_size`) → `save_h5()`.
  - The HDF5 file holds the `data_freq_time` and `data_dm_time` datasets. Its attributes are `cand_id`, `tcand`, `dm`, `dm_opt`, `snr`, `snr_opt`, `width`, `label`, `rfi_mask`, `kill_mask`, `filelist`, `tsamp`, `nchans`, `foff`, `nspectra` and the header fields.

  — [candidate.py](https://github.com/thepetabyteproject/your/blob/main/your/candidate.py), [your_candmaker.py](https://github.com/thepetabyteproject/your/blob/main/bin/your_candmaker.py)
- **RFI** (`utils/rfi.py`). `savgol_filter` fits a Savitzky–Golay filter to the bandpass and flags deviations (Agarwal et al. 2020). `sk_filter` uses spectral kurtosis with median ± σ·1.4826·MAD bounds. `sk_sg_filter` combines the two. — [rfi.py](https://github.com/thepetabyteproject/your/blob/main/your/utils/rfi.py)
- **Other scripts.** `your_heimdall.py` (runs heimdall on a DADA buffer fed by `your` via `utils/heimdall.py:HeimdallManager`), `your_header.py`, `your_viewer.py`, `your_writer.py`, `your_bandpass.py`, `your_rfimask.py`, `your_h5plotter.py`, `your_combine_mocks.py`. — [your bin/](https://github.com/thepetabyteproject/your/tree/main/bin)

#### 10. FETCH (Agarwal et al.; Fast Extragalactic Transient Candidate Hunter)
- **Status and license.** GPL-3.0. The default branch is `tf2` (TensorFlow 2 / Keras). Its last commit is 2025-06-09, although GitHub `pushed_at` is 2026-11-20 (a push to another branch). setup.py says version 0.2.0. — [FETCH repo](https://github.com/devanshkv/fetch)
- **Paper.** Agarwal et al. 2020, MNRAS 497, 1661, "FETCH: A deep-learning based classifier for fast transient classification". — [arXiv:1902.06343](https://arxiv.org/abs/1902.06343)
- **Input.** `your` candidate HDF5 files. The `DataGenerator` (`fetch/data_sequence.py`) reads `data_freq_time` (transposed, detrended, divided by std, median subtracted) and `data_dm_time` (divided by std, median subtracted), each 256×256×1 by default. — [data_sequence.py](https://github.com/devanshkv/fetch/blob/tf2/fetch/data_sequence.py)
- **Models.** These are two-branch networks: one CNN on the frequency–time image (FT) and one on the DM–time image (DMT), with their features merged before a softmax over 2 classes. The 11 pre-trained models a–k are named `<idx>_FT_<net>_<n>_DMT_<net>_<n>_<k>`, for example `a_FT_DenseNet121_2_DMT_Xception_13_256`, `b_FT_DenseNet121_2_DMT_VGG16_1_32`, `c_FT_DenseNet169_6_DMT_Xception_13_112`, … `h_FT_DenseNet201_4_DMT_InceptionResNetV2_20_160`, `k_FT_DenseNet121_2_DMT_InceptionV3_18_64`. Weights are downloaded and verified by MD5 (`get_model`). Fine-tuning goes through `ready_for_train(model, nf, ndt, nft)` with Adam (lr 1e-3) and binary cross-entropy. — [fetch/models](https://github.com/devanshkv/fetch/tree/tf2/fetch/models), [utils.py](https://github.com/devanshkv/fetch/blob/tf2/fetch/utils.py)
- **Output.** `bin/predict.py -m <a..k> -p 0.5` writes `results_<model>.csv` with columns `candidate`, `probability`, `label`. `bin/train.py` handles training. — [predict.py](https://github.com/devanshkv/fetch/blob/tf2/bin/predict.py)

#### 11. ASTROFLOW and DRAFTS (brief)
- **ASTROFLOW** (`lintian233/astroflow`; PyPI `pulseflow`).
  - Status: BSD-3-Clause, last push 2026-06-26.
  - What it is: a real-time end-to-end GPU single-pulse/FRB pipeline. CUDA dedispersion (`src/core/dedispered.cu`) runs in `mode: subband` or the default mode. RFI mitigation uses IQRM (`src/core/iqrm.cu`) and zero-DM (`rfi.use_zero_dm`, `use_iqrm`). A **YOLOv11n** object detector (`modelname: yolov11n`) runs on DM–time images, and there is a boxcar detector (`src/detection/boxcardetection.cpp`) as well as `aidetection.cpp`.
  - Inputs and outputs: it reads `.fil`, PSRFITS and PSRDADA, and writes candidate PNGs and CSVs. Configuration is one YAML file.
  - Paper: "ASTROFLOW: A Real-Time End-to-End Pipeline for Radio Single-Pulse Searches", ApJS, doi:10.3847/1538-4365/ae4a26.

  — [ASTROFLOW repo](https://github.com/lintian233/astroflow), [configuration.md](https://github.com/lintian233/astroflow/blob/main/docs/configuration.md)
- **DRAFTS** (`SukiYume/DRAFTS`).
  - Status: MIT, last commit 2026-10-02.
  - Production stages: (1) FITS organisation and prefetching (FAST-centric; the first two polarisations are summed); (2) GPU dynamic RFI removal, fine-channel dedispersion and time/frequency downsampling (CuPy, with data shared with PyTorch through DLPack); (3) **CenterNet** object detection on 512×512 time–DM images to localise TOA and DM; (4) a **ConvNeXt** binary classifier on dedispersed time–frequency cut-outs.
  - The repo also contains dataset generation, training, and raw8/packed2 injection benchmarks against PRESTO.
  - Paper: arXiv:2410.03200.

  — [DRAFTS repo](https://github.com/SukiYume/DRAFTS)

### Inferences
- The heimdall `.cand` columns (snr, sample, time, filter, dm_idx, dm, members, begin, end) are the de facto interchange format. `your_candmaker`'s `width` column is log₂(samples) precisely so heimdall's `filter_idx` can be passed straight through, and that is how heimdall, your and FETCH became a standard chain.
- heimdall's merging is friends-of-friends with fixed index tolerances (3 time units scaled by width, 3 filter steps, 200 DM trials). TransientX's DBSCAN with `min_samples` additionally rejects isolated noise points. That is the main conceptual difference in duplicate handling.
- riptide's S/N formula and heimdall's 1/√w normalisation are the same matched-filter idea. riptide's version is exact for a folded profile with an unknown baseline, while heimdall assumes baseline removal has already zeroed the mean.

### Gaps
- I did not read the inner FFA kernel (`transforms.hpp`) line by line. The O(m·p·log₂m) description comes from the paper's well-known result and is not quoted from the code.
- I could not locate the pFDMT and TLSM implementation in the PulsarX or XLibs source by name: grep for "fdmt" and "tikhonov" found nothing. The algorithms are confirmed only from the paper abstract, so the file names for them remain unknown.
- The heimdall `clean_filterbank_rfi.cu` statistics (exact zero-DM thresholds per `rfi_tol`) and the dedisp DM-list formula were not transcribed. The dedisp formula is in Barsdell et al. 2012 / the dedisp source, which I did not clone.
- The AA FDAS template generation (z-step math) and the periodicity output format were not read in detail.

## Q2. Where the codes overlap, and where one is clearly best

### Takeaway
For GPU real-time FRB detection, heimdall (mature, simple) and AstroAccelerate (more sensitive SPDT with width-dependent noise statistics, plus FDAS and periodicity in one library) overlap. TransientX is the strongest CPU-only option. For candidate classification, your plus FETCH is the established chain, while ASTROFLOW and DRAFTS fold detection into ML object detectors. For periodicity, peasoup (GPU FFT with time-domain resampling, and now orbital template banks) and riptide (FFA, best at long periods and low duty cycles) complement each other, with PulsarX as the folding and optimisation back end. dmt is the modern, actively developed dedispersion engine (FDMT, direct and coherent; CPU and GPU) that any of these could build on. sigpyproc is the general Python toolbox.

### Cited Findings
- **Dedispersion.**
  - Brute-force GPU (dedisp): heimdall and peasoup. — [heimdall pipeline.cu](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Pipeline/pipeline.cu), [peasoup pipeline_multi.cu](https://github.com/ewanbarr/peasoup/blob/master/src/pipeline_multi.cu)
  - AA's own many-core direct dedispersion (Novotný et al. 2023). — [AA README](https://github.com/AstroAccelerateOrg/astro-accelerate)
  - CPU subband: TransientX. — [arXiv:2401.13834](https://arxiv.org/abs/2401.13834)
  - Pruned FDMT: PulsarX. — [arXiv:2309.02544](https://arxiv.org/abs/2309.02544)
  - FDMT, FDMT-FFT, DDMT, SDMT and CohFDMT on CPU and GPU with exact per-trial σ: dmt. — [dmt README](https://github.com/pravirkr/dmt)
  - Subband dedispersion on GPU: ASTROFLOW. — [configuration.md](https://github.com/lintian233/astroflow/blob/main/docs/configuration.md)
  - Python and numba dedispersion: sigpyproc (`dedisperse`, `tree_dedisperse`). — [base.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/base.py)
- **Matched filtering.**
  - heimdall: power-of-2 boxcars via prefix sum. — [matched_filter.cu](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Pipeline/matched_filter.cu)
  - AA: iterative decimated boxcar bank with width-dependent mean and σ (SPDT; Adámek & Armour 2020). — [aa_device_SPS_long.cu](https://github.com/AstroAccelerateOrg/astro-accelerate/blob/master/src/aa_device_SPS_long.cu)
  - TransientX: geometric widths set by snrloss. — [boxcar.cpp](https://github.com/ypmen/TransientX/blob/master/src/singlepulse/boxcar.cpp)
  - riptide: exact boxcar S/N on folded profiles. — [snr.hpp](https://github.com/v-morello/riptide/blob/master/src/riptide/cpp/snr.hpp)
- **Duplicate removal.**
  - heimdall: GPU friends-of-friends in (time, filter, DM-index). — [are_coincident.cuh](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Pipeline/hd/are_coincident.cuh)
  - TransientX: kd-tree DBSCAN. — [clustering.rst](https://github.com/ypmen/TransientX/blob/master/docs/pipeline/clustering.rst)
  - riptide: frequency clustering plus harmonic flagging. — [pipeline.py](https://github.com/v-morello/riptide/blob/master/src/riptide/pipeline/pipeline.py)
  - peasoup: Harmonic, Acceleration and DM distillers. — [pipeline_multi.cu](https://github.com/ewanbarr/peasoup/blob/master/src/pipeline_multi.cu)
- **RFI.**
  - heimdall: zero-DM narrow and broad excision with random replacement. — [clean_filterbank_rfi.cu](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Pipeline/clean_filterbank_rfi.cu)
  - TransientX and PulsarX: zdot, skewness–kurtosis, kadaneF/T, IQR mask. — [rfi_mitigation.rst](https://github.com/ypmen/TransientX/blob/master/docs/pipeline/rfi_mitigation.rst)
  - sigpyproc: MAD and IQRM. — [core/rfi.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/core/rfi.py)
  - your: Savitzky–Golay and spectral kurtosis. — [rfi.py](https://github.com/thepetabyteproject/your/blob/main/your/utils/rfi.py)
  - AA: zero-DM (with outliers) and rfi flags. — [ska_test_file.txt](https://github.com/AstroAccelerateOrg/astro-accelerate/blob/master/input_files/ska_test_file.txt)
- **Periodicity.**
  - peasoup: time-domain resampling, then FFT and harmonic sum; also Keplerian template banks. — [peasoup README](https://github.com/ewanbarr/peasoup)
  - AA: FDAS with Fourier-domain z (and w) templates, plus a plain periodicity search. — [MANUAL.md](https://github.com/AstroAccelerateOrg/astro-accelerate/blob/master/MANUAL.md)
  - riptide: FFA, best for long periods and narrow duty cycles. — [arXiv:2004.03701](https://arxiv.org/abs/2004.03701)
- **Folding.**
  - PulsarX is the recommended folder for peasoup candidates. — [peasoup README](https://github.com/ewanbarr/peasoup)
  - sigpyproc and riptide also fold. — [base.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/base.py), [time_series.py](https://github.com/v-morello/riptide/blob/master/src/riptide/time_series.py)
- **ML stage.**
  - FETCH classifies FT plus DMT 256×256 images produced by `your_candmaker`. — [data_sequence.py](https://github.com/devanshkv/fetch/blob/tf2/fetch/data_sequence.py)
  - ASTROFLOW runs YOLOv11n on DM–time images. — [configuration.md](https://github.com/lintian233/astroflow/blob/main/docs/configuration.md)
  - DRAFTS runs CenterNet on 512×512 time–DM images, then a ConvNeXt classifier. — [DRAFTS README](https://github.com/SukiYume/DRAFTS)
- **Maintenance snapshot** (Nov 2026):
  - Very active: dmt (v0.8.0, 2026-11-11), PulsarX (2026-11-08), sigpyproc (v2.0.0, 2026-10-17), DRAFTS (2026-10).
  - Active: riptide (2026-09), your (2026-09), TransientX (2026-07), ASTROFLOW (2026-06).
  - Moderate: AA (last commit 2025-12), peasoup (2025-08), FETCH tf2 (2025-06).
  - Low or legacy: heimdall (2023-06).

  — [GitHub API / git log on each repo; heimdall SourceForge](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/)

### Inferences
- Best single-pulse sensitivity per GPU-second is probably AA (its SPDT corrects the noise σ for each width and decimation level), while heimdall is the simplest to deploy and integrate (it is fed by `your_heimdall.py`). For CPU-only clusters, TransientX is clearly the best-documented choice.
- For long-period or narrow pulsars, riptide (FFA) beats FFT codes. For binaries, peasoup's resampling (and template banks) or AA's FDAS and jerk search apply. PulsarX then folds and optimises the candidates.
- For a new Python codebase, sigpyproc (I/O, cleaning, folding and reference implementations of heimdall and PRESTO single-pulse search) plus dmt (dedispersion engines) plus your and FETCH (candidate cut-outs and classification) cover most needs without CUDA C++.

### Gaps
- No head-to-head benchmark across all codes on identical data was found. The comparisons above rest on design and the papers, not on measured sensitivity.

## Q3. Key source files and classes for later reference

### Takeaway
The table below gives the file-level entry points for each code. Repository roots are on GitHub unless noted; heimdall is on SourceForge.

### Cited Findings
- **heimdall.**
  - `Applications/heimdall.C` (main); `Pipeline/pipeline.cu` (`hd_create_pipeline`, `hd_execute`).
  - Pipeline stages: `Pipeline/clean_filterbank_rfi.cu`, `remove_baseline.cu`, `get_rms.cu`, `matched_filter.cu`, `find_giants.cu`, `label_candidate_clusters.cu`, `merge_candidates.cu`, `median_filter.cu`, `measure_bandpass.cu`.
  - Parameters and options: `default_params.C`, `parse_command_line.C`.
  - Headers: `Pipeline/hd/*.h` (`params.h`, `types.h`, `are_coincident.cuh`).
  - Formats: `Formats/SigprocFile.C`, `PSRDadaRingBuffer.C`, `hd/header.h`.
  - Tools: `Applications/coincidencer.C`, `candidate_profiler.C`, `generate_dmlist.C`; `Scripts/trans_gen_overview.py`, `trans_freqplus_plot.py`, `frb_detector.py`.

  — [heimdall tree](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/)
- **AstroAccelerate.**
  - Entry points: `src/aa_main.cpp`; `include/aa_permitted_pipelines_generic.hpp`.
  - API and strategy: `src/aa_py_pipeline_api.cpp`, `aa_ddtr_strategy.cpp`.
  - Dedispersion and binning: `aa_device_dedispersion_kernel.cu`, `aa_corner_turn.cpp`, `aa_bin_gpu.cpp`.
  - Single-pulse search: `aa_device_single_pulse_search*.cu`, `aa_device_SPS_long*.cu`, `aa_device_SPS_inplace*.cu`, `aa_device_MSD*.cu`, `aa_device_analysis.cu`, `aa_device_threshold.cu`, `aa_device_peak_find.cu`.
  - Periodicity and acceleration: `aa_device_periods.cu`, `aa_device_harmonic_summing*.cu`, `aa_fdas_{device,host,util}.cu`, `aa_device_acceleration_fdas.cu`, `aa_device_jerk_search.cu`.
  - Preprocessing: `aa_zero_dm*.cpp`, `aa_device_rfi*.cu`.
  - Input: `aa_sigproc_input.cpp`.

  — [AA src](https://github.com/AstroAccelerateOrg/astro-accelerate/tree/master/src)
- **dmt.**
  - Engines: `lib/fdmt*.cpp`, `lib/ddmt*.cpp`, `lib/sdmt.cpp`, `lib/cfdmt*.cpp`.
  - Search and noise: `lib/boxcar_*.cpp`, `lib/kadane_plan.cpp`, `lib/dm_grids.cpp`, `lib/fdmt_plan_noise.cpp`.
  - Python bindings: `src/dmt_python.cpp`. Public C++ header: `include/dmt/dmt.hpp`.

  — [dmt](https://github.com/pravirkr/dmt)
- **TransientX.**
  - `src/singlepulse/singlepulse.cpp` (`SinglePulse::prepare/run`), `dedisperse_fil.cpp` (main), `boxcar.cpp`, `cluster.cpp`, `candidate.cpp`, `candplot.cpp`, `replot_fil.cpp`, `archivewriter_tx.cpp`.
  - Headers: `include/{dedisperse,boxcar,cluster,clustering,hough_transform,ringbuffer}.h`.
  - Docs: `docs/pipeline/*.rst`.

  — [TransientX](https://github.com/ypmen/TransientX)
- **PulsarX.**
  - `src/pulsar/{filtool,dedisperse_all_fil,psrfold_fil,psrfold_fil2,psrfold_tim,dmffdot,pulsarsearch,gridsearch,dedispersionlite,dedispersionliteU,archivewriter,archivelite,pulsarplot,pulsar_injector,filmaker,psrheader}.cpp`.
  - Fold templates: `include/template/*.template`.
  - Shared library code: the XLibs repo.

  — [PulsarX](https://github.com/ypmen/PulsarX), [XLibs](https://github.com/ypmen/XLibs)
- **riptide.**
  - C++ core: `src/riptide/cpp/{transforms,periodogram,snr,kernels,downsample,running_median,block}.hpp`, bound in `python_bindings.cpp`.
  - Python: `libffa.py`, `search.py` (`ffa_search`), `periodogram.py`, `peak_detection.py`, `clustering.py`, `time_series.py`, `folding.py`, `candidate.py`.
  - Pipeline: `pipeline/{pipeline,dmiter,peak_cluster,worker_pool,config_models}.py`; `reading/{sigproc,presto}.py`.
  - App: `apps/rseek.py`.

  — [riptide](https://github.com/v-morello/riptide)
- **peasoup.**
  - Pipelines: `src/pipeline_multi.cu` (Worker, DMDispenser), `pipeline.cpp`, `pipeline_single.cpp`.
  - Kernels: `kernels.cu` (resampling, interbinning, harmonic sums), `folding_kernels.cu`.
  - Transforms: `include/transforms/*.hpp`. Data types: `include/data_types/*.hpp`.
  - Utilities: `include/utils/utils.hpp` (AccelerationPlan), `output_stats.hpp`, `xml_util.hpp`.

  — [peasoup](https://github.com/ewanbarr/peasoup)
- **sigpyproc.**
  - Core classes: `src/sigpyproc/{base,readers,header,timeseries,fourierseries,foldedcube,block}.py`.
  - Kernels: `kernels/*.py` (numba).
  - Modules: `core/{rfi,stats,dispersion,filters,candidates}.py`, `algorithms/{box,sps,heimdall,period}.py`, `io/{sigproc,pfits,fbh5,digifil}.py`.
  - CLI: `apps/spp_*.py`.

  — [sigpyproc3](https://github.com/FRBs/sigpyproc3)
- **your.**
  - Classes: `your/your.py` (Your, Header), `your/writer.py` (Writer), `your/candidate.py` (Candidate).
  - Format back-ends: `your/formats/*.py`. Utilities: `your/utils/{rfi,heimdall,gpu,plotter}.py`.
  - Scripts: `bin/your_*.py`.

  — [your](https://github.com/thepetabyteproject/your)
- **FETCH.**
  - `fetch/{data_sequence,utils}.py`, `fetch/models/<a..k>_FT_*_DMT_*/*.json` plus `model_list.csv`, and `bin/{predict,train}.py`.

  — [FETCH](https://github.com/devanshkv/fetch/tree/tf2)
- **ASTROFLOW and DRAFTS.**
  - ASTROFLOW: `src/core/{dedispered.cu,iqrm.cu}`, `src/detection/{aidetection,boxcardetection}.cpp`, `src/rfi/rfimarker.cu`, `src/data/{filterbank,psrfits,psrdada}.cpp`, `python/astroflow/{search,frbdetector,dmtime,cli}.py`.
  - DRAFTS: `search_pipeline/`, `object_detection/`, `binary_classification/`, `dataset_generation/`, `injection_benchmark/`.

  — [ASTROFLOW](https://github.com/lintian233/astroflow), [DRAFTS](https://github.com/SukiYume/DRAFTS)

### Inferences
- For a flowchart-level port, the minimal file set is:
  - Single-pulse: heimdall `pipeline.cu`, `matched_filter.cu`, `remove_baseline.cu`, `are_coincident.cuh` and `label_candidate_clusters.cu`; TransientX `singlepulse.cpp`, `boxcar.cpp` and `cluster.cpp`.
  - Periodicity: peasoup `pipeline_multi.cu`, `kernels.cu` and `utils.hpp`; riptide `periodogram.hpp` and `snr.hpp`.

### Gaps
- Exact file locations of PulsarX's pFDMT and TLSM (they are probably in XLibs under different names) were not identified.
- The ASTROFLOW and DRAFTS internals were surveyed only at README and configuration level.
