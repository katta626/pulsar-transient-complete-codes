# DSPSR, PSRCHIVE and pulsar timing software (TEMPO, TEMPO2, PINT): pipeline, formats and math

Method note: most findings come from reading the current source trees, using shallow clones taken in Nov 2026:
- DSPSR `git.code.sf.net/p/dspsr/code`, HEAD 3f3d5b4, 2026-05-27, `AC_INIT([DSPSR],[2026-03-31])`
- PSRCHIVE `git.code.sf.net/p/psrchive/code`, HEAD d77bc07, 2026-04-17, `AC_INIT([PSRCHIVE],[2026-02-26])`
- TEMPO2 `bitbucket.org/psrsoft/tempo2`, HEAD b74248f, 2026-04-21
- PINT `github.com/nanograv/PINT`, HEAD d88d1d2, 2026-11-17

The TEMPO polyco spec was fetched from the TEMPO reference manual. Several doc pages returned 404 or redirected: the dspsr.sourceforge.net manual, the psrchive pat manual and the RTD architecture page. In those cases the source code (usage strings, `set_help` text) is cited instead.

Source-link conventions:
- DSPSR: `https://sourceforge.net/p/dspsr/code/ci/master/tree/<path>`
- PSRCHIVE: `https://sourceforge.net/p/psrchive/code/ci/master/tree/<path>`
- TEMPO2: `https://bitbucket.org/psrsoft/tempo2/src/master/<path>`
- PINT: `https://github.com/nanograv/PINT/blob/master/<path>`

---

## DSPSR: operation chain, formats, coherent dedispersion, filterbank, detection, SK, folding, CUDA, CLI

### Takeaway
DSPSR is a C++ library of `dsp::Operation` / `dsp::Transformation<In,Out>` classes. The core data path is:

File/Input → `BitSeries` → Unpacker → `TimeSeries` → [Filterbank / Convolution with a `Dedispersion` response] → `Detection` → `Fold` → `PhaseSeries` → `Archiver`, which writes a PSRCHIVE `Pulsar::Archive`.

Coherent dedispersion multiplies the spectrum by the cold-plasma chirp. The chirp phase is computed in `Dedispersion::build_compute`, and the transform is done with overlap-save. The FFT length is the power of two ≥ 2× the smearing length (in samples) that minimises nfft·log nfft/(nfft−nfilt). Folding uses PSRCHIVE `Pulsar::Predictor` objects, either a TEMPO polyco or a TEMPO2 2-D Chebyshev predictor, generated from a par file (`-E`) or loaded directly (`-P`).

### Cited Findings

**Identity, licence and citation**
- DSPSR is "an open-source C++ library for digital signal processing of pulsar timeseries". It provides "coherent dedispersion, filterbank formation, pulse folding", reads about 40 formats including PSRFITS, VDIF and Mark5, and has docs that moved to ReadTheDocs — [dspsr.readthedocs.io](https://dspsr.readthedocs.io/en/latest/); [dspsr.sourceforge.net/manuals/dspsr](http://dspsr.sourceforge.net/manuals/dspsr/)
- The licence is the Academic Free License v2.1 (`LICENSE`). `CITATION` asks users to cite van Straten & Bailes (2011), http://www.publish.csiro.au/?paper=AS10021 (PASA 28, 1) — [DSPSR LICENSE / CITATION](https://sourceforge.net/p/dspsr/code/ci/master/tree/CITATION)

**Source layout**
- `Kernel/Classes`: the base data containers and I/O, with these files:
  - `Operation.C`, `Input.C`, `File.C`, `MultiFile.C`, `Unpacker.C`
  - `BitSeries.C`, `TimeSeries.C`, `WeightedTimeSeries.C`, `DataSeries.C`, `Observation.C`
  - `TwoBitCorrection.C` (2-bit dynamic-level setting, used for "2-bit excision" via `-2`)
  - `FourBitUnpacker.C`, `EightBitUnpacker.C`, `FloatUnpacker.C`, `GenericEightBitUnpackerCUDA.cu`
  - `MemoryCUDA.C`, `TransferCUDA.C`, `TimeSeriesCUDA.cu`, `DADAFile.C`, `DADAHeader.C`, `ASCIIObservation.C`
  - [Kernel/Classes](https://sourceforge.net/p/dspsr/code/ci/master/tree/Kernel/Classes/)
- `Kernel/Formats`: one directory per backend or format.
  - apsr, asp, bcpm, boreas, bpsr, caspsr, cpsr, cpsr2
  - dada, dummy, emerlin, fadc, fits (`FITSFile.C`, `FITSUnpacker.C`, `GUPPIFITSUnpacker.C`, `FITSOutputFile.C`, `FITSDigitizer.C`)
  - gmrt, guppi, kat, lbadr, lbadr64, lofar_dal, lump, lwa
  - mark4, mark5, mark5b, maxim, mopsr, mwa, pdev, pmdaq, puma, puma2
  - s2, sigproc (`SigProcFile.C`, `SigProcObservation.C`, `SigProcUnpacker.C`, `SigProcOutputFile.C`, `SigProcDigitizer.C`)
  - ska1, spda1k, spigot, uwb, vdif, wapp
  - Registered in `File_registry.C` and `Unpacker_registry.C` — [Kernel/Formats](https://sourceforge.net/p/dspsr/code/ci/master/tree/Kernel/Formats/)
- `Signal/General`: the DSP operations.
  - Convolution and dedispersion: `Convolution.C` (+`ConvolutionCUDA.cu`, `ConvolutionCUDASpectral.cu`), `Response.C`, `PlasmaResponse.C`, `Dedispersion.C`, `DedispersionPipe.C`
  - Filterbanks: `Filterbank.C` (+`FilterbankEngineCPU.C`, `FilterbankEngineCUDA.cu`), `InverseFilterbank.C` (+CUDA), `PolyPhaseFilterbank.C`, `TFPFilterbank.C`, `IncoherentFilterbank.C`
  - Detection and spectral kurtosis: `Detection.C` (+`DetectionCUDA.cu`), `SpectralKurtosis.C`, `SKFilterbank.C`, `SKDetector.C`, `SKMasker.C` (all with CUDA versions)
  - Other operations: `Apodization.C`, `OptimalFFT.C`, `Rescale.C`, `TScrunch.C`, `FScrunch.C`, `PScrunch.C`, `Derotation.C`, `PolnCalibration.C`, `SampleDelay.C`
  - Pipeline classes: `SingleThread.C`, `MultiThread.C`, `LoadToFil.C`, `LoadToFITS.C`, `LoadToQuantize.C`
  - Companion tools: `digifil.C`, `digifits.C`, `digihist.C`, `digistat.C`, `dmsmear.C`
  - [Signal/General](https://sourceforge.net/p/dspsr/code/ci/master/tree/Signal/General/)
- `Signal/Pulsar`: `Fold.C` (+`FoldCUDA.cu`), `PhaseSeries.C`, `Archiver.C`, `ArchiverExtensions.C`, `LoadToFold1.C`, `LoadToFoldN.C`, `LoadToFoldConfig.C`, `TimeDivide.C`, `Subint.h`, `CyclicFold.C` (+CUDA), `PhaseLockedFilterbank.C`, `FoldManager.C`, `UnloaderShare.C`, `t2pred.C`, and the `dspsr.C` main program — [Signal/Pulsar](https://sourceforge.net/p/dspsr/code/ci/master/tree/Signal/Pulsar/)

**Class wiring**
- `Fold` is declared as `Transformation <TimeSeries, PhaseSeries> ("Fold", outofplace)` (Fold.C line 27). Its output is a `PhaseSeries`, which `Archiver` converts to a `Pulsar::Archive` — [Signal/Pulsar/Fold.C](https://sourceforge.net/p/dspsr/code/ci/master/tree/Signal/Pulsar/Fold.C)
- `LoadToFold::Config::add_options` chains three option groups: `SingleThread::Config::add_options`, the `dedisp` (DedispersionPipe) options, and the alias "add dspsr -t nthread alias", i.e. `menu.find("threads")->set_short_name('t')` — [Signal/Pulsar/LoadToFoldConfig.C](https://sourceforge.net/p/dspsr/code/ci/master/tree/Signal/Pulsar/LoadToFoldConfig.C)

**Coherent dedispersion math**
- The dispersion constant is coded as `dm_dispersion = 2.41e-4`. The source comment is `DM (pc cm^-3) = 2.410331(2)×10^-4 D (s MHz^2)`.
- Delay: `delay_time(freq) = (DM/2.41e-4) / freq^2` seconds (freq in MHz).
- Transfer-function phase, from `build_compute(chan_freq f0, freq f)` (f is the offset from channel centre):
  - `coeff = -sign(bw) · 2π · (1e6·DM/2.41e-4) / f0^2`
  - `phase(f) = coeff · f^2/(f0+f) + delay_phase`
  - `delay_phase = -2π f · interchannel_fractional_sample_delay` corrects the fractional-sample inter-channel delay.
  - Phasors are built as `polar(1.0, phase)`, i.e. H(f) = exp(i·phase).
  - With `build_delays` set, phase/(2πf) gives delays in µs instead.
  - [Signal/General/Dedispersion.C](https://sourceforge.net/p/dspsr/code/ci/master/tree/Signal/General/Dedispersion.C)
- Smearing:
  - `smearing_time(cfreq,bw) = delay_time(cfreq-|bw|/2) - delay_time(cfreq+|bw|/2)`.
  - The smear length is computed per half of the band, giving `impulse_pos` and `impulse_neg`.
  - It is padded by `smearing_buffer = 0.1` ("add another ten percent… pollution due to the cyclical convolution effect"), or by 5% in psrdisp-compatible mode.
  - `nsmear = ceil(tsmear·sampling_rate)`. A threshold of 16·1024·1024/nchan samples forces a symmetric impulse response.
  - `impulse_pos` samples are dropped from the beginning of each cyclical-convolution result and `impulse_neg` from the end. This is the overlap-save discard.
  - [Signal/General/PlasmaResponse.C](https://sourceforge.net/p/dspsr/code/ci/master/tree/Signal/General/PlasmaResponse.C)
- FFT length choice (`OptimalFFT::get_nfft(nfilt)`):
  - Start at `nfft_min` = the smallest power of 2 strictly greater than `nfilt`.
  - Iterate powers of 2 up to a benchmark maximum.
  - Minimise cost T(nfft) = nfft·log2(nfft)/(nfft−nfilt), since "ngood = nfft − nfilt".
  - `-x nfft|minX` overrides this. `--fft-bench` uses measured benchmarks.
  - [Signal/General/OptimalFFT.C](https://sourceforge.net/p/dspsr/code/ci/master/tree/Signal/General/OptimalFFT.C)

**Dispersion-removal options** (`DedispersionPipeConfig.C`)
- `-F <N>[:D]` "create an N-channel filterbank". The `:D` suffix requests coherent dedispersion inside the filterbank (the synthetic filterbank).
- `-IF <N>[:D]` "create inverse filterbank with N output channels".
- `-dr` deripple.
- `-t-taper` / `-f-taper` set the apodization window.
- `-D dm` "dispersion measure".
- `-K` "remove inter-channel dispersion delays".
- `-x nfft|minX`.
- `-R` "apply time-variable narrow-band RFI filter".
- `-pac dbase` "calibrator database for phase-coherent matrix convolution".
- `-derotate` / `-rm` for phase-coherent Faraday rotation correction.
- `-fft-bench`.
- [Signal/General/DedispersionPipeConfig.C](https://sourceforge.net/p/dspsr/code/ci/master/tree/Signal/General/DedispersionPipeConfig.C)

**Detection**
- `-d npol` help text: "0=P,Q; 1=PP+QQ; 2=PP,QQ; 3=(PP+QQ)^2; 4=PP,QQ,PQ,QP".
- `Detection.C` states: `Intensity`, `PPQQ`, `NthPower`, `PP_State`, `QQ_State`, `Coherence` ("PP, QQ, Re[PQ], Im[PQ]"), `Stokes` ("Stokes I,Q,U,V").
- `stokes_detect.ic` computes:
  - S0 = pp+qq
  - S1 = pp−qq
  - S2 = 2(p_r q_r + p_i q_i)
  - S3 = 2(p_r q_i − p_i q_r)
  - where pp = |p|², qq = |q|² for the two complex polarisations p, q.
- `-4` computes fourth-order moments. `-n ndim` is "[experimental] ndim of output when npol=4".
- [Signal/General/Detection.C](https://sourceforge.net/p/dspsr/code/ci/master/tree/Signal/General/Detection.C); [stokes_detect.ic](https://sourceforge.net/p/dspsr/code/ci/master/tree/Signal/General/stokes_detect.ic); [dspsr.C](https://sourceforge.net/p/dspsr/code/ci/master/tree/Signal/Pulsar/dspsr.C)

**Spectral kurtosis**
- The SK estimator in `SKFilterbank.C` is `SK = ((M+1)/(M-1)) · (M·S2/S1² − 1)`, where S1 = Σ P_k and S2 = Σ P_k² over M power samples per channel and polarisation.
- Time-scrunched and frequency-scrunched SK are also formed.
- `SpectralKurtosis.C` sets the default `std_devs(3.0)`.
- [Signal/General/SKFilterbank.C](https://sourceforge.net/p/dspsr/code/ci/master/tree/Signal/General/SKFilterbank.C); [SpectralKurtosis.C](https://sourceforge.net/p/dspsr/code/ci/master/tree/Signal/General/SpectralKurtosis.C)
- SK CLI options:
  - `-skz` "apply spectral kurtosis (SK) based RFI mitigation"
  - `-skzm samples`, `-skzover`, `-skzs stddevs`, `-skz_start/-skz_end chan`
  - `-skz_config sk.yaml`, `-skz_no_fscr`, `-skz_no_tscr`, `-skz_no_ft` (despeckler)
  - `-skzn threads`, `-sk_fold`, `-noskz_too` ("also produce un-zapped version of output")
  - [dspsr.C](https://sourceforge.net/p/dspsr/code/ci/master/tree/Signal/Pulsar/dspsr.C)

**Folding and predictors**
- `Fold::get_folding_predictor(params, observation)` works as follows:
  - Calls `Pulsar::Generator::factory(params)`. For a `Tempo::Predict` (polyco) generator it sets `nspan`/`ncoef`.
  - Sets site = telescope, time span = start−864 s to +86400 s (comment: "make a predictor valid for the next 24 hours… a bit of a hack"), and frequency span = centre ± bw/2.
  - Calls `generate()`.
- `Fold::maximum_nbin = 1024`, `minimum_bin_width = 1.2` input samples, `power_of_two = true`. `choose_nbin()` picks nbin from period and time resolution unless `-b` is given.
- When only a predictor is present, period = 1/`folding_predictor->frequency(start_time)`.
- [Signal/Pulsar/Fold.C](https://sourceforge.net/p/dspsr/code/ci/master/tree/Signal/Pulsar/Fold.C)

**Key dspsr options** (`set_help` strings in `dspsr.C` / `SingleThread.C`)
- `-E file` "pulsar ephemeris used to generate predictor"
- `-P file` "phase predictor used for folding"
- `-w file` "phase predictors used for folding."
- `-b nbin` "number of phase bins in folded profile"
- `-c period` "folding period (in seconds)", with `-cepoch MJD` and `-p phase`
- `-X name` "additional pulsar to be folded"
- `-asynch-fold` "fold on CPU while processing on GPU"
- `-A` "output single archive with multiple integrations"
- `-nsub N`
- `-s` "create single pulse sub-integrations"
- `-turns N`
- `-L seconds` "create integrations of specified duration", with `-Lepoch`, `-Lmin`, `-y`
- `-e ext`, `-O name`, `-fname`, `-Y` (pdmp extras)
- `-j job` / `-J a.psh` (psrsh commands run before output)
- `-G nbin` phase-locked filterbank; `-cyclic N` / `-cyclicoversample M` cyclic spectra
- `-B`, `-f`, `-k`, `-N` override bandwidth, frequency, telescope and source name
- `-C offset`, `-m MJD`
- `-scloffs` "denormalize PSRFITS using DAT_SCL and DAT_OFFS"
- `-2 code` "2-bit excision"
- `-S seek`, `-T total`, `-seek MJD`
- `-set key=value`, `-list`, `-header`, `-W` (disable weights), `-backends`
- `-minram`, `-U MB|minX` "upper limit on RAM usage"
- `-threads nthread` (alias `-t`), `-cpu cores`
- `-cuda devices` "comma-separated list of CUDA devices"
- `-Z lib` FFT library
- `-r` timing report
- `-dump op`, `-mask op`
- Sources: [dspsr.C](https://sourceforge.net/p/dspsr/code/ci/master/tree/Signal/Pulsar/dspsr.C); [SingleThread.C](https://sourceforge.net/p/dspsr/code/ci/master/tree/Signal/General/SingleThread.C)

**CUDA, Python and companion tools**
- CUDA is built in through `config/cuda.m4` and `config/Makefile.cuda`. CUDA engines exist for Convolution, Filterbank, InverseFilterbank, Detection, Fold, CyclicFold, SK (Computer, Detector, Filterbank, Masker), Rescale, the scrunches, SampleDelay, and unpackers (Float, GenericEightBit, FITSDigitizer) — [Signal/General](https://sourceforge.net/p/dspsr/code/ci/master/tree/Signal/General/); [config](https://sourceforge.net/p/dspsr/code/ci/master/tree/config/)
- A SWIG Python binding exists: `python/dspsr.i` and `example.py` — [dspsr/python](https://sourceforge.net/p/dspsr/code/ci/master/tree/python/)
- `t2pred.C` is a standalone TEMPO2-predictor test or utility in `Signal/Pulsar` — [t2pred.C](https://sourceforge.net/p/dspsr/code/ci/master/tree/Signal/Pulsar/t2pred.C)
- Companion programs share this pipeline:
  - `digifil` writes SIGPROC filterbanks (`LoadToFil`)
  - `digifits` writes PSRFITS search mode (`LoadToFITS`)
  - `dmsmear` computes smearing
  - [Signal/General](https://sourceforge.net/p/dspsr/code/ci/master/tree/Signal/General/)

### Inferences
- The transfer function in the code matches the standard form H(f0+f) = exp[∓i 2π D·DM f²/(f0²(f0+f))] with D = 1/2.41×10⁻⁴ ≈ 4149 s MHz² pc⁻¹ cm³ (Hankins & Rickett 1975). The sign is set by the sideband (`sign(bw)`). Frequencies are in MHz, and the factor 1e6 converts to Hz·s units.
- The overlap-save "nfilt" is impulse_pos + impulse_neg (the dispersive smear in samples × 1.1). The minimum FFT is the next power of two above nfilt, and the efficiency kept is (nfft−nfilt)/nfft.
- A synthetic filterbank (`-F N:D`) does coherent dedispersion in each output channel. A large forward FFT is followed by N small backward FFTs of length nfft/N, and the per-channel chirp is applied in between. This reduces the FFT length needed per channel, compared with a single full-band convolution followed by a filterbank.
- The SK estimator form matches Nita & Gary (2010). Its expected value is 1 for Gaussian noise, with variance ≈ 4/M for large M. The ±`skzs`σ (default 3) thresholds are applied around 1.

### Gaps
- The van Straten & Bailes (2011) paper itself could not be retrieved: the arXiv IDs tried were wrong and the PASA page was not fetched. Paper-specific equations and figures (e.g. the exact synthetic-filterbank description) are inferred from code, not quoted from the paper.
- The dspsr.readthedocs.io application-manual and architecture subpages returned 404 at their guessed URLs, so official prose descriptions of options were not obtained. The `set_help` strings were used instead.
- The exact list of the "~40 formats" was taken from the `Kernel/Formats` directory names. Which of them are enabled by default at configure time was not checked.

---

## PSRCHIVE: class hierarchy, plugins, tools, TOA algorithms, calibration, Python, licence

### Takeaway
PSRCHIVE's base library (`Base/Classes`) defines `Pulsar::Archive` (abstract) → `Integration` (sub-integration) → `Profile` (one polarisation × channel). Concrete file formats are plugins under `Base/Formats`, loaded through `Archive::Agent` / `Archive::Advocate<Plugin>`.

`pat` measures TOAs by template matching, choosing among PGS (default Fourier phase gradient, Taylor 1992), GIS, PIS, ZPS, SIS, FDM, COF and RVM, plus matrix template matching (`-p`). Output is in parkes, tempo2, itoa or princeton format. The licence is AFL 2.1. Python bindings are SWIG plus the NumPy C-API. No explicit NumPy-2 statement was found.

### Cited Findings

**Identity, licence and citation**
- `CITATION` asks users to cite Hotan, van Straten & Manchester (2004), http://www.publish.csiro.au/?paper=AS04022 (PASA 21, 302).
- `COPYING` is the Academic Free License v2.1.
- The latest configure version is `2026-02-26`.
- [PSRCHIVE CITATION](https://sourceforge.net/p/psrchive/code/ci/master/tree/CITATION)

**Class hierarchy and format plugins**
- `Base/Classes` contains `Archive.C`, `Integration.C`, `Profile.C`, `BasicArchive.C`, `BasicIntegration.C`, `Agent.C`, `Archive_load.C`, `Archive_unload.C`, `IntegrationManager.C`, `ProfileAmps.C`, `MoreProfiles.C`, `FourthMoments.C`, `Calibrator.C`, `FITSAlias.C` and others.
- `Archive.h` declares "Advocates the use of an Archive derived class: `class Agent;`" and "Advocates the use of Plugin: `template<class Plugin> class Advocate;`".
- [Base/Classes](https://sourceforge.net/p/psrchive/code/ci/master/tree/Base/Classes/)
- Format plugins in `Base/Formats`:
  - ASCII, ASP, BPP, EPN, EPOS
  - PRESTO
  - PSRFITS, including the `psrheader.fits` template that defines the PSRFITS fold-mode header
  - PuMa
  - Timer (the original Parkes/CPSR `.ar`/`.rf`/`.cf` format)
  - UVM, WAPP, pdv_KtZ
  - [Base/Formats](https://sourceforge.net/p/psrchive/code/ci/master/tree/Base/Formats/)

**Applications** (`More/Applications`)
- Archive manipulation: `pam`, `paz`, `pazi`, `psradd`, `psrsplit`, `psrsh`, `psrzap`, `psrwt`, `psrtrash`
- Inspection and plotting: `pdv`, `psrplot`, `pav`, `psrstat`, `psrtxt2`, `print_Predictor`
- Timing: `pat`
- Polarimetric calibration: `pac`, `pacv`, `pcm` (+`pcmavg`, `pcmdiff`, `pcmerr`, `pcmrm`, `pcmtxt`)
- Flux: `fluxcal`, `psrflux`
- Rotation measure and polarisation: `rmfit`, `psrpol`, `psr4th`
- Analysis: `paas`, `pas`, `pdmp`, `psrspa`, `psrpca`, `psrsmooth`, `psrwhite`, `psrdiff`, `psrmodel`, `psrstability`, `smint`, `spa`, `mtm`, `zap_periodic_spikes`, `dynamic_spectra`, `drifting_subpulses`
- [More/Applications](https://sourceforge.net/p/psrchive/code/ci/master/tree/More/Applications/)

**pat options** (usage text in `pat.C`)
- Preprocessing: `-F` fscrunch, `-T` tscrunch, `-j/-J` jobs.
- Templates and harmonics:
  - `-s stdfile` standard; `-a stdfiles` auto-select a standard
  - `-n harmonics`; `-c` choose the maximum harmonic
  - `-P` "Do not fscrunch the standard (frequency-dependent template)"
  - `-I nsubband` "Output nsubband TOAs with delta-DM"
  - `-D` denoise the standard
- `-S period` zaps periodic spikes (implies SIS).
- `-p` "Enable matrix template matching"; `-o` outputs the MTM model parameters.
- `-A name`, "default: PGS". The available algorithms are:
  - "PGS = Fourier phase gradient"
  - "GIS = Gaussian interpolation"
  - "PIS = Parabolic interpolation"
  - "ZPS = Zero pad interpolation"
  - "SIS = Sinc interpolation of cross-corration function"
  - "FDM = Fourier domain with Markov chain Monte Carlo"
  - "COF = Centre of Flux"
  - "RVM = Rotating Vector Model magnetic meridian"
- `-m` paas component model; `-g` Gaussian model.
- `-f "fmt <flags>"` "[default: parkes] Available formats: parkes tempo2, itoa, princeton". The tempo2 flags are `i`, `r` and `IPTA`.
- `-X "text"` adds extra TOA-line text. `-C`, `-r`, `-R`, `-u`, `-b` and `-K/-t` control other output and plots.
- [More/Applications/pat.C](https://sourceforge.net/p/psrchive/code/ci/master/tree/More/Applications/pat.C)
- The pat manual itself is at https://psrchive.sourceforge.net/manuals/pat. Per a search snippet, it says the PGS method is used if no algorithm is specified with -A — [psrchive pat manual](https://psrchive.sourceforge.net/manuals/pat)

**TOA estimator source files** (`More/Timing`)
- `ArrivalTime.C`, `ArrivalTime_tempo2.C`
- `PhaseGradShift.C` (PGS), `FourierDomainFit.C` (FDM), `GaussianShift.C` (GIS), `ParIntShift.C` (PIS), `ZeroPadShift.C` (ZPS), `SincInterpShift.C` (SIS), `FluxCentroid.C` (COF), `RotatingVectorModelShift.C` (RVM)
- `MatrixTemplateMatching.C` (MTM), `ComponentModel.C`, `MeanArrivalTime.C`, `ShiftEstimator.C`
- [More/Timing](https://sourceforge.net/p/psrchive/code/ci/master/tree/More/Timing/)
- `ScalarTemplateMatching.h` documents the core method: "This class implements the Taylor (1992) algorithm, using the Van Wijngaarden–Dekker–Brent method to minimize the objective merit function". It was originally the function `model_profile(npts, narrays, prf, std, scale, sigma_scale, shift, sigma_shift, chisq, …)` (still in `More/General/model_profile.C`) — [ScalarTemplateMatching.h](https://sourceforge.net/p/psrchive/code/ci/master/tree/More/General/Pulsar/ScalarTemplateMatching.h)
- Copyright headers credit Willem van Straten (e.g. `MatrixTemplateMatching.C` 2009–2017) — [MatrixTemplateMatching.C](https://sourceforge.net/p/psrchive/code/ci/master/tree/More/Timing/MatrixTemplateMatching.C)
- `pcm` documents a related but separate "Measurement Equation Template Matching (METM)" method (van Straten 2013) used for calibration modelling — [pcm manual METM](https://psrchive.sourceforge.net/manuals/pcm/metm.shtml)
- A comparative study of TOA methods (PGS, FDM and others) is Wang et al. 2022, A&A 658, A181 (arXiv:2111.13482) — [arXiv:2111.13482](https://arxiv.org/abs/2111.13482)

**Calibration classes** (`More/Polarimetry`)
- Flux: `FluxCalibrator*.C` (ConstantGain, VariableGain, Policy, Ext) and `FluxCalManager.C`
- Polarisation: `HybridCalibrator.C`, `ManualPolnCalibrator.C`, plus the SingleAxis / ReceptionModel classes used by `pac`/`pcm`
- [More/Polarimetry](https://sourceforge.net/p/psrchive/code/ci/master/tree/More/Polarimetry/)

**Python bindings**
- They are SWIG (`More/python/psrchive.i`).
- They include `numpy/arrayobject.h` and call `import_array()`.
- They return NumPy arrays via `PyArray_SimpleNewFromData` / `PyArray_SimpleNew(…, NPY_FLOAT)`, e.g. `get_amps`, a weighted-amplitude copy, and a 4-D data copy.
- A `PSRCHIVE_PyArray_DATA` macro wraps `PyArray_DATA`.
- The fallback `setup.py` still uses `from distutils.core import setup`.
- [More/python/psrchive.i](https://sourceforge.net/p/psrchive/code/ci/master/tree/More/python/psrchive.i)
- NumPy 2.0 breaks C-API binary compatibility, so extensions built against NumPy 1.x will not load and must be recompiled against 2.x — [NumPy 2.0 migration guide](https://numpy.org/doc/2.0/numpy_2_0_migration_guide.html)

### Inferences
- PGS math (Taylor 1992):
  - With profile P and template S DFTs P_k, S_k, minimise χ²(a,τ) = Σ_{k=1}^{n_h} |P_k − a S_k e^{−2πikτ}|²/σ².
  - Equivalently, maximise Σ |P_k||S_k| cos(φ_k − θ_k + 2πkτ). Here a is the scale factor, τ the phase shift in turns, and σ the per-harmonic noise.
  - The uncertainty comes from ∂²χ²/∂τ².
  - The `model_profile` signature (scale, shift, their sigmas, chisq) matches this.
- MTM generalises the scalar fit to the full Stokes (or coherency) profiles, with a Jones-matrix transformation. This is consistent with van Straten (2006, 2013) and with the `-o` "best-fit matrix model parameters" option.
- The pat default output is the "parkes" format. For TEMPO2/PINT use, `-f tempo2` (optionally with `IPTA`) gives `FORMAT 1` lines with flags.
- The PSRCHIVE Python bindings use only long-stable NumPy C-API calls. They should compile against NumPy 2 when rebuilt, but a wheel or module built against NumPy 1.x will fail to import under NumPy 2. The distutils-based `setup.py` fallback will not work on Python ≥3.12, so the autotools build path must be used.

### Gaps
- No primary source gave a PSRCHIVE-specific NumPy-2 status statement (release note, issue or commit). The shallow clone had no history for `More/python`.
- The Hotan et al. 2004 and van Straten, Demorest & Oslowski 2012 (AR&T 9, 237) papers were not fetched. Statements about them rely on the CITATION file and general knowledge.
- The pat manual page body could not be fetched (404 on the direct URL). Exact PGS equations are therefore reconstructed from Taylor (1992) knowledge rather than quoted.
- Option-level details for pam, paz, pazi, psradd, pdv, psrplot, psrstat, pac, pcm, psrflux and rmfit were not extracted. Only their existence and source files were confirmed.

---

## Timing software: TEMPO, TEMPO2, PINT — models, par/tim, barycentring, polyco and predictor math, fitting

### Takeaway
Polycos (TEMPO format, also produced by TEMPO2 `-polyco` and PINT `pint.polycos`) give phase = RPHASE + 60·DT·F0 + Σ COEFF(n)·DT^(n−1), with DT in minutes from TMID.

TEMPO2 predictors (`t2pred.dat`, produced by `tempo2 -pred`) are piecewise 2-D Chebyshev series in (MJD, frequency) plus an explicit DISPERSION_CONSTANT/ν² term. They are evaluated by `ChebyModel_GetPhase` in `cheby2d.c`, and both DSPSR and PSRCHIVE consume them via `Pulsar::Predictor`.

PINT (BSD licence, Python ≥3.9, numpy ≥1.23) provides WLS/GLS/Downhill/LM/Wideband fitters. Its latest release found was 1.1.8 (2026-10-08).

### Cited Findings

**TEMPO polyco format**
- Line 1:
  - pulsar name (cols 1–10)
  - date dd-mmm-yy (11–19)
  - UTC hhmmss.ss (20–31)
  - TMID MJD (32–51)
  - DM (52–72)
  - Doppler ×10⁻⁴ (74–79)
  - log10 fit RMS in periods (80–86)
- Line 2:
  - RPHASE (1–20)
  - F0 reference rotation frequency (21–38)
  - observatory number (39–43)
  - span in minutes (44–49)
  - NCOEFF (50–54)
  - observing frequency in MHz (55–75)
  - binary phase (76–80)
- Then coefficient lines of 3 × 25 columns.
- Formulas:
  - DT = (T − TMID)·1440
  - PHASE = RPHASE + DT·60·F0 + COEFF(1) + DT·COEFF(2) + DT²·COEFF(3) + …
  - FREQ(Hz) = F0 + (1/60)(COEFF(2) + 2·DT·COEFF(3) + 3·DT²·COEFF(4) + …)
- [TEMPO polyco spec](https://tempo.sourceforge.net/ref_man_sections/tz-polyco.txt)
- PINT implements the same formulas:
  - `evalabsphase` computes `dt = (t - tmid)*MIN_PER_DAY` and `phase += rphase + Phase(dt*60.0*f0)`.
  - Frequency is `f0 + s/60.0`.
  - The docstrings state ΔT = 1440(T−T_mid).
  - [src/pint/polycos.py](https://github.com/nanograv/PINT/blob/master/src/pint/polycos.py)

**TEMPO2 command line**
- `tempo2 -polyco "mjd1 mjd2 nspan ncoeff maxha sitename freq"` writes polycos.
- `-pred "sitename mjd1 mjd2 freq1 freq2 ntimecoeff nfreqcoeff seg_length (s)"` "Creates a predictive 2D Chebyshev polynomial".
- Environment variables: `$TEMPO2` (T2runtime: clocks, observatories, EOPs), `$TEMPO2_CLOCK_DIR`, `$TEMPO2_PLUG_PATH`.
- [tempo2 help.txt](https://bitbucket.org/psrsoft/tempo2/src/master/help.txt)
- In polyco mode, if TZRMJD is unset, TEMPO2 warns "tzrmjd not set. Setting to" PEPOCH — [tempo2.C](https://bitbucket.org/psrsoft/tempo2/src/master/tempo2.C)

**TEMPO2 predictor (`t2pred.dat`)**
- File structure: each segment is written as a block:
  - `ChebyModel BEGIN`
  - `PSRNAME`
  - `SITENAME`
  - `TIME_RANGE mjd_start mjd_end`
  - `FREQ_RANGE f_start f_end`
  - `DISPERSION_CONSTANT`
  - `NCOEFF_TIME nx`
  - `NCOEFF_FREQ ny`
  - `COEFFS` lines, three per line
  - `ChebyModel END`
- Phase evaluation:
  - x = −1 + 2(mjd−mjd_start)/(mjd_end−mjd_start)
  - y = −1 + 2(freq−f_start)/(f_end−f_start)
  - phase = Cheby2D_Evaluate(coeffs, x, y) + dispersion_constant/freq²
- Frequency comes from a derivative series, scaled by 2/(mjd_end−mjd_start)/86400 to give Hz.
- The basis is T_ij(x,y) = T_i(x)T_j(y). Coefficients are c_ij = Σ_a Σ_b f(x_a,y_b) T_ij(x_a,y_b) at the Chebyshev nodes x = cos(π(k−½)/n), and evaluation uses a 2-D Clenshaw recurrence with a ½·c_0 term.
- `tempo2pred.c` dispatches `T2Predictor_GetPhase` / `GetFrequency` to either `ChebyModelSet` or `T1PolycoSet`. TEMPO1 polycos are a second predictor kind in the same API, and `T2Predictor_Insert` is "not implemented for T1 polyco".
- [cheby2d.c](https://bitbucket.org/psrsoft/tempo2/src/master/cheby2d.c); [tempo2pred.c](https://bitbucket.org/psrsoft/tempo2/src/master/tempo2pred.c)

**TEMPO2 source layout and licence**
- The licence is GPL v3 (`COPYING`).
- Key files:
  - par/tim reading: `readParfile.C`, `readTimfile.C`
  - clock and time scales: `clkcorr.C`, `toa2utc.C`, `utc2tai.C`, `tai2tt.C`, `tt2tdb.C`
  - barycentring and delays: `formBats.C`, `get_obsCoord.C`, `shapiro_delay.C`, `dm_delays.C`, `sw_delay.C`, `tropo.C`, `jpleph.c`, `readEphemeris*.C`
  - residuals and fitting: `formResiduals.C`, `doFit.C`, `t2fit*.C`, `cholesky*.C`
  - polyco and predictors: `polyco.C`, `t1polyco.c`
  - binary models: `BTmodel.C`, `DDmodel.C`, `DDKmodel.C`, `DDGRmodel.C`, `DDSmodel.C`, `DDHmodel.C`, `ELL1model.C`, `ELL1Hmodel.C`, `ELL1kmodel.C`, `MSSmodel.C`, `T2model.C`
  - [tempo2 repo](https://bitbucket.org/psrsoft/tempo2/src/master/)

**PINT**
- Timing-model components (`src/pint/models`):
  - `astrometry.py`, `spindown.py`, `dispersion_model.py`, `solar_wind_dispersion.py`, `solar_system_shapiro.py`, `troposphere_delay.py`
  - `binary_bt.py`, `binary_dd.py`, `binary_ddk.py`, `binary_ell1.py`
  - `glitch.py`, `jump.py`, `phase_offset.py`, `absolute_phase.py` (TZR), `wave.py`, `wavex.py`, `dmwavex.py`, `cmwavex.py`, `chromatic_model.py`, `fdjump.py`, `ifunc.py`, `piecewise.py`
  - `noise_model.py`, `timing_model.py`, `model_builder.py`
  - [PINT models](https://github.com/nanograv/PINT/tree/master/src/pint/models)
- Fitters (`fitter.py`): `Fitter`, `DownhillFitter`, `DownhillWLSFitter`, `DownhillGLSFitter`, `WidebandDownhillFitter`, `PowellFitter`, `WLSFitter`, `GLSFitter`, `WidebandTOAFitter`, `LMFitter`, `WidebandLMFitter` — [fitter.py](https://github.com/nanograv/PINT/blob/master/src/pint/fitter.py)
- `get_TOAs(…, ephem, include_bipm, bipm_version, planets, …)` sets barycentring options. Planet Shapiro terms use ("jupiter","saturn","venus","uranus","neptune","earth"), and solar-system positions come from `objPosVel_wrt_SSB` — [toa.py](https://github.com/nanograv/PINT/blob/master/src/pint/toa.py)
- PINT licence and requirements: BSD, `requires-python >=3.9`, `numpy>=1.23.0`, `astropy>=5.0.5` — [pyproject.toml](https://github.com/nanograv/PINT/blob/master/pyproject.toml)
- CHANGELOG releases: 1.1.8 (2026-10-08), 1.1.7 (2026-08-24), 1.1.6 (2026-08-10), 1.1.5 (2026-05-06), 1.1.4 (2025-08-18), 1.1.3 (2025-04-25) — [CHANGELOG.md](https://github.com/nanograv/PINT/blob/master/CHANGELOG.md)
- PINT v1.0.0 was archived on Zenodo on 2024-04-26 — [Zenodo 11396241](https://zenodo.org/records/11396241)
- The accompanying paper is Susobhanan et al. 2024, "PINT: Maximum-likelihood estimation of pulsar timing noise parameters", ApJ 971, 150 (arXiv:2405.01977) — [arXiv:2405.01977](https://arxiv.org/abs/2405.01977); [INSPIRE record](https://inspirehep.net/files/976a8d8127e68bdfece7627ba20f8088)

### Inferences
- **Timing-model chain**, common to all three packages:
  1. Site arrival time in UTC(obs).
  2. Clock corrections: observatory clock → GPS → UTC → TAI → TT(BIPM).
  3. Time-scale conversion: TT → TDB (or TCB in TEMPO2 by default).
  4. Barycentring: Roemer delay −r_obs·n̂/c (from the JPL DE ephemeris), solar Shapiro delay, and Einstein delay. The TT→TDB step corresponds to the Einstein delay.
  5. Dispersion delay D·DM/ν² with the TEMPO2/PINT convention D = 1/2.41×10⁻⁴ s MHz², plus solar-wind and troposphere terms.
  6. Binary delays.
  7. Spin phase φ(t) = φ0 + F0·Δt + F1·Δt²/2 + F2·Δt³/6 + … at PEPOCH.
  8. Residual = φ − nearest integer, divided by F0.
  9. Linearised least squares (WLS, or GLS with a noise covariance) solves for parameter offsets from the design matrix ∂φ/∂p.
  - Hobbs, Edwards & Manchester (2006, MNRAS 369, 655) and Edwards et al. (2006, MNRAS 372, 1549) document this for TEMPO2. Those papers were not fetched in this session.
- **par/tim conventions** (from general knowledge, not fetched):
  - A par file holds lines of `PARAM value [fitflag] [uncertainty]`: PSRJ, RAJ, DECJ, F0, F1, PEPOCH, DM, BINARY, TZRMJD, TZRFRQ, TZRSITE, EPHEM, CLK, UNITS.
  - TEMPO2 tim files start with `FORMAT 1`, then lines `file freq(MHz) MJD unc(µs) site -flag value …`. This matches pat's `-f tempo2` output.
- **How the files flow:**
  - `.par` → `dspsr -E` → `Pulsar::Generator` produces a polyco or t2pred predictor. The predictor is stored in the output PSRFITS fold-mode archive (POLYCO or T2PREDICT binary table) together with the ephemeris (PSRPARAM table).
  - `pat -s std` then turns the `.ar` into `.tim` TOAs, and tempo2 or PINT fits the `.par` to the `.tim`.

### Gaps
- The TEMPO2 papers (Hobbs et al. 2006; Edwards et al. 2006), Luo et al. 2021 (ApJ 911, 45, PINT) and the ATNF TEMPO2 documentation pages were not fetched. Statements on barycentring components and par/tim syntax are therefore general knowledge, flagged as inference.
- Which PSRFITS binary tables (POLYCO, T2PREDICT, PSRPARAM) PSRCHIVE's PSRFITS plugin reads and writes was not verified in source (`Base/Formats/PSRFITS`) in this session.
- Exact TEMPO (TEMPO1) source structure and licence were not examined. Only its polyco spec was fetched.
