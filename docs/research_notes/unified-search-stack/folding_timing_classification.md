# Folding, Archive Handling, Timing, and Candidate Classification: Best Existing Implementations

Scope: which existing code to adopt, wrap, or reimplement for a unified stack (Rust CPU kernels, C++/CUDA GPU kernels, Python orchestration). Research date: late Nov 2026. Note on method: web access was limited (several primary pages returned HTTP 403, e.g. A&A full text and psrchive.sourceforge.net manuals). Anything not tied to a fetched or searched source is placed under Inferences or Gaps and labelled as prior knowledge.

## Which parts of DSPSR and PSRCHIVE are irreplaceable vs which can be reimplemented faster? What is their GPU support status (dspsr CUDA)?

### Takeaway
For search-candidate folding, PulsarX (psrfold_fil) is the strongest existing implementation: about 50x faster than PRESTO prepfold on large candidate sets, with FDMT-based dedispersion and a regularised least-squares folding algorithm. It is the natural model for a Rust or CUDA folding kernel. DSPSR remains the reference for coherent dedispersion and folding of baseband with full polarimetry, and it has real multi-GPU CUDA support. PSRCHIVE's calibration, polarimetry, and TOA (pat) machinery is the part that is hardest to replace. Simple archive I/O and manipulation (pam/pdv/psradd/paz-style scrunching, zapping, and summing) is straightforward to reimplement.

### Cited Findings
**Folding: PulsarX**
- PulsarX (Men et al. 2023, A&A 679, A20; arXiv:2309.02544) is a complete folding pipeline for large-scale surveys, covering RFI mitigation, dedispersion, folding, and parameter optimisation. — [ADS abstract](https://ui.adsabs.harvard.edu/abs/2023arXiv230902544M/abstract); [arXiv PDF](https://arxiv.org/pdf/2309.02544)
- It uses a pruned variant of the fast discrete dispersion measure transform (FDMT) and a new folding algorithm based on Tikhonov-regularised least squares. — [ADS abstract](https://ui.adsabs.harvard.edu/abs/2023arXiv230902544M/abstract)
- Benchmark: psrfold_fil2 is about 50x faster than prepfold with many candidates. The paper's benchmark names `psrfold_fil2`, while survey papers name `psrfold_fil`, so the exact binary name should be checked in the repo. — [arXiv 2309.02544](https://arxiv.org/pdf/2309.02544)
- TRAPUM (MeerKAT) SMC survey uses PulsarX psrfold_fil to cut the compute cost of large candidate numbers. — [TRAPUM SMC survey, arXiv:2405.12029](https://arxiv.org/pdf/2405.12029)
- PulsarX is listed in the ASCL. — [ASCL/ADS](https://ui.adsabs.harvard.edu/abs/2023ascl.soft12012M/abstract)

**DSPSR (baseband, coherent dedispersion, GPU)**
- dspsr is C++ and CUDA, described as a high-performance, general-purpose tool for high-time-resolution radio pulsar studies using accelerated computing. CHIME/Pulsar uses it in production. — [CHIME Pulsar system overview, arXiv:2008.05681](https://arxiv.org/pdf/2008.05681)
- CHIME/Pulsar runs coherent dedispersion and timing-model folding on GPU buffers before writing to disk. dspsr runs much faster on (8+8)-bit input. — [arXiv:2008.05681](https://arxiv.org/pdf/2008.05681)
- The dspsr paper (van Straten & Bailes 2011) lists phase-coherent dispersion removal with optional Jones-matrix convolution for polarimetry, and CUDA extension. Fold output is PSRCHIVE-readable archives. — [dspsr paper (PASA)](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/187B3FC44A3610F092A337FF6A7CB0A6/S1323358000000606a.pdf/dspsr-digital-signal-processing-software-for-pulsar-astronomy.pdf); [arXiv:1008.3973](https://arxiv.org/pdf/1008.3973)
- CUDA build: set PACKAGES to the CUDA install and run configure. Usage: `-cuda 0,0,1,1` starts one POSIX thread per entry, and repeated IDs add CUDA streams on that device to overlap PCIe transfers with compute. Docs have moved to ReadTheDocs (dspsr.readthedocs.io). The old manual page refers to CUDA 4.2-era examples. — [dspsr GPU manual](https://dspsr.sourceforge.net/manuals/dspsr/gpu.shtml)
- The SKA PST (Pulsar Timing) team has an MR "added dspsr cuda support" in ska-pst-dsp-tools, so dspsr CUDA is being built into SKA-era tooling. — [SKA GitLab MR](https://gitlab.com/ska-telescope/pst/ska-pst-dsp-tools/-/merge_requests/6)
- Other GPU real-time coherent dedispersion systems exist (GMRT pipeline; GUPPI/NRAO GPU benchmarking). They are useful as references for kernel design. — [GMRT real-time coherent dedispersion, arXiv:1509.00186](https://arxiv.org/pdf/1509.00186); [NRAO GPU benchmarking](https://www.cv.nrao.edu/~pdemores/gpu)
- dspsr project activity is on SourceForge. — [DSPSR SourceForge activity](https://sourceforge.net/p/dspsr/activity)

**PSRCHIVE**
- PSRCHIVE is an open-source C++ library for analysing pulsar data. It is the format and analysis backbone that dspsr depends on. — [PSRCHIVE project](https://psrchive.sourceforge.net/); [dspsr paper](https://arxiv.org/pdf/1008.3973)
- Timing practice still relies on PSRCHIVE `pat`. Its Fourier-domain MCMC (FDM) and Matrix Template Matching (MTM) methods appear in recent PTA and timing papers, often combined with PulsePortraiture templates. — [PulsePortraiture GitHub](https://www.github.com/pennucci/PulsePortraiture); [J1528-3146 timing, arXiv:2304.06578](https://arxiv.org/pdf/2304.06578); [PPTA DR3, arXiv:2306.16230](https://arxiv.org/pdf/2306.16230)

### Inferences
- Treat these as irreplaceable, so wrap rather than rewrite:
  - PSRCHIVE's polarimetric calibration (Pulsar/Measurement-Equation modelling, pac/pcm) and the pat TOA estimators (FDM/MTM), which PTA results depend on.
  - dspsr's baseband unpackers for many backends, plus its coherent-dedispersion-with-Jones-matrix path.
  These encode decades of instrument-specific correctness that would cost a lot to re-validate.
- Treat these as reimplementable faster, in Rust or CUDA:
  - Search-mode folding of filterbank candidates. Copy the PulsarX design: batch many candidates, subband, FDMT, regularised fold.
  - Archive scrunching and zapping (pam/paz equivalents), summing (psradd), text dumping (pdv), and plotting (psrplot).
  - Reading and writing PSRFITS fold-mode as plain FITS binary tables.
- For a search pipeline (rather than baseband timing), DSPSR is not needed on the hot path. PulsarX-style folding of filterbank plus a PSRFITS writer covers candidate folding. Keep dspsr as an optional backend for baseband or coherent follow-up.
- The fetched sources do not show GPU folding inside PulsarX, which appears to be CPU-optimised. That leaves room for a CUDA batch-folding kernel, and for a Metal kernel if Apple Silicon is a target.

### Gaps
- I could not confirm which dspsr operations are GPU-accelerated. From prior knowledge (unverified this session): the filterbank, convolution (coherent dedispersion), detection, folding, and some unpackers are, and SpectralKurtosis has a CUDA path. The fetched GPU manual page says nothing on this, and ReadTheDocs was not fetched.
- I found no dspsr 2025–2026 release notes. CUDA 12.x compatibility is not confirmed from a primary source.
- I could not access GPU folding in peasoup or sigpyproc's fold implementation details this session. From prior knowledge (unverified): peasoup does GPU dedispersion and FFT acceleration search but hands folding to other tools, and sigpyproc's fold is a simple CPU (numba/C) fold over filterbank.
- No head-to-head benchmark of PulsarX vs dspsr vs sigpyproc folding beyond the 50x prepfold figure.

## Maintenance status 2025–2026, build pain points, licenses (GPL/AFL), Python bindings

### Takeaway
PSRCHIVE is actively maintained: a release dated 2026-02-26, commits in late 2025, and conda-forge packaging including aarch64. However, its Python bindings are documented as incompatible with NumPy 2.x, which is a major friction point for a modern Python orchestration layer. FETCH is GPL-3.0 on TensorFlow 2. PINT is the maintained Python-native timing package (v1.0 paper in 2024). Licences for PSRCHIVE and dspsr were not confirmed from a primary source.

### Cited Findings
- PSRCHIVE's SourceForge activity shows a release dated 2026-02-26, a commit on 25 Nov 2025, and conda-forge work including "Add aarch64 on conda-forge". The git remote is git.code.sf.net/p/psrchive/code. — [PSRCHIVE SourceForge activity](https://sourceforge.net/p/psrchive/activity)
- conda-forge hosts psrchive builds. — [anaconda.org conda-forge/psrchive](https://anaconda.org/conda-forge/psrchive/files/manage?page=2)
- The PSRCHIVE Python interface ships with the standard distribution. Per search-result snippets of the manual, the bindings are incompatible with NumPy 2.x. The manual page itself returned 403, so this is from the search snippet only. — [PSRCHIVE Python manual](https://psrchive.sourceforge.net/manuals/python)
- Pure-Python alternatives exist:
  - PyPulse (mtlam) is a pure-Python package for handling and analysing PSRFITS. — [PyPulse GitHub](https://github.com/mtlam/PyPulse); [PyPI](https://pypi.org/project/PyPulse/0.1.0/)
  - PSRpy (Fonseca) is another option. — [PSRpy GitHub](https://github.com/emmanuelfonseca/PSRpy)
- PulsePortraiture requires PSRCHIVE (Python bindings), lmfit, numpy, and scipy, so it inherits the PSRCHIVE/NumPy constraint. — [PulsePortraiture GitHub](https://www.github.com/pennucci/PulsePortraiture)
- PINT v1.0.0 was described in a 2024 paper that adds maximum-likelihood single-pulsar noise estimation with AIC model comparison. — [PINT noise paper, arXiv:2405.01977](https://arxiv.org/html/2405.01977v2)
- The original PINT paper (Luo et al. 2021) validated PINT against Tempo 13.101 and Tempo2 2019.01.1 on PSR J1600-3053 (12,433 TOAs, NANOGrav 11-yr). — [PINT paper, arXiv:2012.00074](https://arxiv.org/pdf/2012.00074)
- FETCH (devanshkv/fetch, tf2 branch) is GPL-3.0 and on TensorFlow ≥2, with 71 commits, no GitHub releases, and a Zenodo citation at v0.1.8. Candidates come as HDF5 made by `your_candmaker.py` from the `your` package. — [FETCH GitHub](https://github.com/devanshkv/fetch)

### Inferences
- Python orchestration should avoid importing PSRCHIVE's SWIG bindings in the main environment because of the NumPy 2 conflict. Options:
  1. A native PSRFITS reader/writer (astropy.io.fits or fitsio, or a Rust crate with PyO3), with PyPulse as a reference.
  2. Running PSRCHIVE CLI tools (pat, pac, pam) in an isolated conda env or container as subprocesses.
- Licensing: FETCH is GPL-3.0. Vendoring or linking GPL code (and PSRCHIVE/dspsr if GPL) into a permissively licensed unified codebase would impose GPL terms. Calling them as separate processes or optional plugins keeps the core licence clean.
- From prior knowledge (unverified this session): PSRCHIVE and dspsr are understood to be under the Academic Free License (AFL) v2.1 or v3.0. Tempo2 and PRESTO are GPL. PINT and sigpyproc are BSD/MIT-style permissive. Verify each repo's COPYING/LICENSE before deciding.
- Build pain points (prior knowledge, unverified):
  - PSRCHIVE and dspsr use autotools with many optional dependencies: FFTW, CFITSIO, tempo2, PGPLOT, and on the dspsr side psrdada and CUDA.
  - dspsr must be configured against a matching PSRCHIVE build, and backend selection goes through a `backends.list` file.
  - conda-forge packaging (now with aarch64) is the least painful route, including on Apple Silicon. CUDA does not apply on macOS.

### Gaps
- PSRCHIVE/dspsr licence text was not confirmed from a primary source. The search did not surface the COPYING file.
- I could not verify whether newer PSRCHIVE (2026 release) fixed NumPy 2 compatibility.
- No primary data found on the latest PINT/tempo2 release versions for 2025–2026, or on tempo2's maintenance status.
- Ephemeris and predictor handling was not deeply verified this session. From prior knowledge (unverified): tempo1 uses polycos, tempo2 uses Chebyshev 2-D phase predictors (T2predict), and PINT can generate polycos and has a native Python phase model. dspsr and PSRCHIVE call tempo or tempo2 to make predictors, which is a runtime dependency to plan for.

## Best-performing open classifiers with published metrics, and whether they run on CPU/GPU/Apple Silicon

### Takeaway
Published metrics are almost all self-reported on the authors' own test sets, and independent checks show sharp drops. For pulsar candidates:
- PICS (2014) is the historical baseline and has been retrained for FAST.
- Newer CNN/DenseNet models (2025) report ms-scale inference.

For single-pulse/FRB candidates:
- FETCH is still the de facto open standard (>99.5% self-reported accuracy/recall).
- An independent 2025 study found much lower recall on injected signals.
- DRAFTS (2024/2025, open on GitHub) reports >99% accuracy and finds more than 3x the bursts Heimdall found on FAST data.
- A 2026 transformer preprint claims to beat DRAFTS and RaSPDAM.

All are standard PyTorch/TensorFlow CNNs, so they run on CPU or GPU, and on Apple Silicon via TF-Metal or PyTorch MPS. I found no published Apple-Silicon benchmarks.

### Cited Findings
**Single-pulse / FRB**
- FETCH (Agarwal et al. 2020) provides 11 models, each with accuracy and recall above 99.5% on a test set of real RFI and pulsar candidates. — [FETCH arXiv:1902.06343](https://arxiv.org/abs/1902.06343v2)
- In the Northern HTRU single-pulse paper (2025), overall recall from FETCH on injected signals was much lower than the >0.99 reported, for all models. The authors note no accuracy was deduced, so the recall is only indicative. — [N-HTRU II, arXiv:2511.17797](https://arxiv.org/pdf/2511.17797) (snippet from search result; abstract fetch did not show the numbers)
- DRAFTS (Zhang et al.) is an object-detection stage plus a binary classifier. On FRB 20190520B FAST data it found more than 3x the bursts Heimdall found. Its GitHub claims >99% classification accuracy on FAST and GBT data. — [DRAFTS arXiv:2410.03200](https://arxiv.org/html/2410.03200v1); [DRAFTS GitHub](https://github.com/SukiYume/DRAFTS)
- Connor & van Leeuwen (2018) was an early deep-learning single-pulse classifier using dynamic spectra, DM-time, and multi-beam inputs. — [arXiv:1803.03084](https://arxiv.org/pdf/1803.03084)
- MeerTRAP uses an ML classifier (FRBID) for MeerKAT transients. — [MeerTRAP FRB sample, arXiv:2302.10107](https://arxiv.org/pdf/2302.10107)
- Newer work:
  - A "DM-free" ML FRB search on multibeam data. — [arXiv:2512.19249](https://arxiv.org/pdf/2512.19249)
  - A 2026 transformer-based framework reported to outperform PRESTO/Heimdall and the RaSPDAM/DRAFTS baselines on accuracy and speed. — [arXiv:2603.05958](https://arxiv.org/pdf/2603.05958) (per search summary; not fetched)
  - Frabjous (2025), an FRB morphology classifier, not an RFI/real classifier. — [arXiv:2507.14854](https://web3.arxiv.org/pdf/2507.14854) (per search summary)

**Pulsar candidates (folded)**
- PICS (Zhu et al. 2014) is an image-pattern-recognition classifier using deep neural nets, trained on PALFA candidates. — [PICS arXiv:1309.0776](https://arxiv.org/abs/1309.0776)
- PICS has been retrained and customised for the FAST drift-scan survey, with labelled FAST data published. — [arXiv:1903.06383](https://web3.arxiv.org/abs/1903.06383)
- A semi-supervised GAN classifier (2020) used a retrained PICS as its baseline. — [arXiv:2010.07457](https://arxiv.org/pdf/2010.07457)
- Fu et al. 2025 sift "snapshot" candidates from folded dedispersed time series with deep learning. They report about a 10x reduction in folding time and 100% recall on NGC 5904 benchmark data, fine-tuned on MWA and tested on GBT. — [arXiv:2511.06711](https://arxiv.org/html/2511.06711v1)
- A multi-scale DenseNet hybrid (2025) has about 2.3M parameters and 4.2 ms inference per candidate. — [MDPI Appl. Sci. 15(24) 13097](https://www.mdpi.com/2076-3417/15/24/13097)
- A survey-agnostic CNN uses four diagnostic plots with GAN augmentation (thesis). — [arXiv:2510.25774](https://arxiv.org/pdf/2510.25774)

### Inferences
- Recommended architecture: keep the classifier as a pluggable Python stage with ONNX export, so that CPU, CUDA, and Apple Silicon (CoreML/MPS) all run from one model file.
  - FETCH's GPL-3.0 licence and TF2 dependency argue for retraining an equivalent small CNN (two inputs: frequency-time and DM-time) in PyTorch, rather than vendoring.
  - FETCH-compatible HDF5 candidate files (via `your`) remain useful as a shared interchange format.
- Pulsar classifiers:
  - PICS is old (Python 2 era, Theano-like stack per prior knowledge) and the search found no Python 3 release.
  - A modern stack would train a compact CNN on PulsarX/prepfold diagnostic planes (profile, phase-time, phase-frequency, DM curve), PICS-style.
  - Survey-specific retraining is consistently needed, as with the FAST PICS work and the MWA→GBT fine-tuning.
- The independent N-HTRU result shows self-reported >99% figures cannot be trusted for low-S/N regimes. The unified stack should include an injection-recovery harness to measure classifier recall in-situ.

### Gaps
- I found no primary sources this session for SPINN (Morello et al. 2014), the Lyon/GH-VFDT classifiers (Lyon et al. 2016), the GHRSS classifier, SPEGID, CHIME/FRB's ML RFI classifier, or the TRAPUM-specific ML classifier. Their metrics, licences, and code availability are unverified here.
- From prior knowledge (unverified): TRAPUM/MeerKAT pipelines have used PICS-derived models and PulsarX's candidate output, and PulsarX ships a classifier helper. This needs confirmation.
- No published Apple-Silicon (MPS/CoreML) benchmarks for any of these classifiers were found.
- Exact per-model FETCH recall numbers from arXiv:2511.17797 were not retrieved (the abstract lacks them).
- DRAFTS licence and framework (presumed PyTorch) not confirmed.
