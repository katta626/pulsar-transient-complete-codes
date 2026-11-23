# Periodicity (pulsar) search implementations: FFT/accel/jerk, FFA, template-bank, sifting (state as of Nov 2026)

Scope: choosing a "best existing implementation" per function (dereddening, FFT, Fourier interpolation/interbinning, acceleration/jerk search, harmonic summing, time-domain resampling, FFA, template-bank binary search, sifting) for a unified stack with Rust CPU kernels, C++/CUDA GPU kernels and Python orchestration.

Repo metadata (license / last push / stars) below was pulled from the GitHub REST API on 2026-11-24 (`gh api repos/<owner>/<repo>`); each is cited to the repo URL.

## Published sensitivity and speed comparisons

### Takeaway
On sensitivity, the order is now: coherent harmonic summing (Ransom 2026, CoherentSearch.jl / PRESTO `coherentsearch`) > FFA (riptide) > incoherent FFT+harmonic-sum (accelsearch). For speed on binaries, GPU jerk/accel search (PrestoZL, AstroAccelerate) runs about 30–90x faster than CPU PRESTO. But CPU PRESTO `accelsearch` itself got about 8x faster in late 2026, so old GPU-vs-CPU speedup figures overstate the gap.

### Cited Findings
**FFA vs FFT (Morello et al. 2020, riptide)**
- Morello et al. find the FFA "closely approaches the theoretical optimum sensitivity to all periodic signals" and is significantly more sensitive than FFT + incoherent harmonic summing (IHS), whatever the period and duty cycle. They recommend that surveys include an FFA search, and say riptide has found sources that existing FFT+IHS codes missed — [arXiv:2004.03701](https://arxiv.org/abs/2004.03701)
- The paper adds a "search efficiency" term to the radiometer equation. The average FFA efficiency is E_FFA ≈ 0.930 for Gaussian pulses searched with boxcar filters whose trial widths are spaced by factors of 2. The FFT+IHS efficiency is shown only as a function of duty cycle and number of harmonics (h = 1–32) in Fig. 3, with no single number. The FFA's advantage grows at shorter duty cycles — [ar5iv 2004.03701](https://ar5iv.labs.arxiv.org/html/2004.03701)
- Interbinning has a worst-case S/N loss of 7.4% (FFT efficiency divided by 1.074). A double-peaked pulse matched by single-boxcar filters can drop to about 71% (1/√2) — [ar5iv 2004.03701](https://ar5iv.labs.arxiv.org/html/2004.03701)
- riptide speed: one 9-minute DM trial takes about 1 s on a 16-core node, with ≥1000 phase bins for P ≥ 1 s. The benchmark was 128 DM trials of 2^23 samples at 64 µs, searching 1–120 s on Xeon Gold 6140 (OzStar). Scaling with core count was sublinear (memory-bandwidth bound). Complexity is about T·ν0²·b³·log(T·ν0·b), so cost becomes large for P ≲ 500 ms — [ar5iv 2004.03701](https://ar5iv.labs.arxiv.org/html/2004.03701)
- The riptide paper does not cover FFA acceleration searches; the method is for isolated or long-period pulsars — [ar5iv 2004.03701](https://ar5iv.labs.arxiv.org/html/2004.03701)

**Coherent harmonic summing (new in 2026)**
- Ransom (arXiv:2610.11749, Oct 2026), "A Coherent Harmonic Summing Pulsar Search Code" (Julia, CoherentSearch.jl):
  - It rebuilds each trial frequency's profile from about 60 interpolated complex harmonics and tests boxcar templates on it.
  - On more than 860,000 injected pulsars at matched false-alarm rate in white noise, it recovered 76%, against 64% for riptide (recommended configuration) and 42% for PRESTO accelsearch.
  - Its threshold does not shift in red noise.
  - It runs 1.4–2.2x faster than riptide on one core at matched frequency coverage, up to 27x faster on 48 cores, and on modern GPUs up to 9x faster than a 20-core workstation.

  [arXiv:2610.11749](https://arxiv.org/abs/2610.11749)
- PRESTO v6.2.0 adds `coherentsearch`, a C port of CoherentSearch.jl, plus `coherent_sift.py`:
  - It interpolates each harmonic's complex amplitude exactly, inverse-FFTs the harmonic stack into a profile and scores it with boxcar matched filters.
  - `-maxdecim` does harmonic decimation, and a `-deep` pass uses 120 harmonics.
  - It has no acceleration search.
  - Its output is byte-identical to the Julia code on the test file.

  [PRESTO CHANGELOG](https://github.com/scottransom/presto/blob/master/CHANGELOG.md)

**Jerk search**
- Andersen & Ransom 2018 (ApJL, doi:10.3847/2041-8213/aad59f) extend Fourier-domain acceleration search to a linearly changing acceleration (jerk) inside PRESTO. It improves sensitivity when the observation spans about 5–15% of the orbital period, at "significantly longer runtime". It found Ter5am in archival GBT data — [arXiv:1807.07900](https://arxiv.org/abs/1807.07900)
- PRESTO's notes: jerk search "makes searches take a *lot* longer". Typically set `-wmax` to 3–5x `-zmax`, and `-zmax` rarely needs to exceed 300 — [PRESTO CHANGELOG v2.1](https://github.com/scottransom/presto/blob/master/CHANGELOG.md)

**GPU accelsearch / jerk speedups**
- PrestoZL (Zhejiang Lab; Mao et al. 2025, ApJS 280, 36): GPU jerk search in `accelsearch_cu`. On a 1816.77 s M5 observation with zmax=200 and wmax=300 it took 9115 min with PRESTO C, 178 min with PrestoZL (51.19x) and 161 min with the pipelined version (56.38x), on one A40 and a 20-core CPU. Results are claimed to be "fully identical" to PRESTO v4.0, including detections, parameters and SNR. Its GPU `prepsubband_cu` is said to be about 100x faster than CPU `prepsubband` — [PrestoZL README](https://github.com/zhejianglab/PrestoZL)
- A 2026 FAST paper reports 56.4x end-to-end speedup for accelsearch_cu on an A100, compared with OpenMP CPU accelsearch — [arXiv:2601.01912](https://arxiv.org/pdf/2601.01912)
- AstroAccelerate GPU jerk search, on a Titan V, is "almost 90x" faster than 4-thread parallel PRESTO. Harmonic summing was excluded from the timing, though they note it costs at least as much as the matched filtering. Mean/stdev computation and candidate selection dominated their implementation — [arXiv:1911.01353](https://arxiv.org/pdf/1911.01353)
- Laidler's (UCT) GPU accelsearch PhD thesis reports 30–70x faster than serial CPU PRESTO, and says it allowed SKA1-mid survey data to be searched in real time — [UCT thesis](https://open.uct.ac.za/bitstreams/0f91b6df-7927-400f-9aeb-af56390dda85/download)
- GPU FDAS (Dimoudi et al.) uses the correlation technique in AstroAccelerate — [arXiv:1804.05335](https://export.arxiv.org/pdf/1804.05335); improved version at [arXiv:1711.10855](https://arxiv.org/pdf/1711.10855)
- GPU harmonic sum (Adámek & Armour) benchmarks several harmonic-sum variants on CPU and GPU, with the greedy variant about 6.4x faster than a standard GPU baseline — [arXiv:1812.02647](https://arxiv.org/pdf/1812.02647). Adámek, Roy & Armour later describe a greedy 1-D harmonic sum algorithm with a sensitivity evaluation — [arXiv:2202.12817](https://arxiv.org/pdf/2202.12817)
- An INAF 2025 technical talk on porting PRESTO jerk search to GPU says jerk-enabled searches "can take up to a week" — [INAF Perugia 2025 slides](https://indico.ict.inaf.it/event/3173/contributions/23112/attachments/10332/21248/Perugia%20technical%20meeting%202025.pdf)

**CPU PRESTO accelsearch is now much faster (this shifts the GPU baselines)**
- PRESTO v6.3.3 (Nov 2026), with `accelsearch -zmax 200 -numharm 8` on a 100M-bin spectrum on a 16-core Sapphire Rapids:
  - On one thread it is about 1.25x faster than v6.3.2, 1.6x faster than v5.0.0, and about 8x faster than v5.3.1–v6.1.1. The v5.3.1–v6.1.1 slowdown came from six out-of-line `fmaf()` calls in the inner loop on non-FMA builds.
  - On 8 threads it is 7.7x faster than v5.0.0.
  - A `-march=native` build gives about 8% more.
  - Harmonic-sum stages now add all subharmonics with the threshold scan row by row while the row is in cache (`add_search_ffdotpows`), and the subharmonic gather is vectorised.
  - Rednoised input is now searched without per-block median normalisation. On noise this gives 6–11% fewer false candidates above sigma 2.

  [PRESTO CHANGELOG](https://github.com/scottransom/presto/blob/master/CHANGELOG.md)
- In v6.1.2, PRESTO's own benchmark jerk search (`-zmax 50 -wmax 200`) dropped from 297 s to 126 s — [PRESTO CHANGELOG](https://github.com/scottransom/presto/blob/master/CHANGELOG.md)

**peasoup throughput**
- In the HTRU-S medlat/hilat GPU reprocessing (peasoup-type pipeline), median times per DM trial were 43 ms for dedispersion, 128 ms for acceleration search and 182 ms for the full search — [arXiv:1811.04929](https://arxiv.org/pdf/1811.04929)
- TRAPUM/MeerKAT runs peasoup (time-domain resampling + FFT + IHS) on every coherent and incoherent beam, on APSUSE: 60 nodes with 2x GTX 1080 Ti each, up to 864 beams, about 46 TiB/hour for 480 beams. Candidates are folded with PulsarX — [arXiv:2212.08528](https://arxiv.org/pdf/2212.08528), [arXiv:2602.20811](https://arxiv.org/pdf/2602.20811)
- In MMGPS, PulsarX `psrfold_fil` cut folding time "from months to days" — [URSI AP-RASC 2025](https://www.ursi.org/proceedings/procAP25/papers/0943.pdf); PulsarX paper [arXiv:2309.02544](https://arxiv.org/pdf/2309.02544)

**Template-bank / Einstein@Home**
- Balakrishnan et al. (arXiv:2112.11991) describe a coherent GPU search over 3 (circular) or 5 (elliptical) Keplerian parameters using random/stochastic template banks, following Messenger 2009, Harry 2009, Knispel 2011 and Allen 2013 — [arXiv:2112.11991](https://arxiv.org/pdf/2112.11991)
- The Einstein@Home radio search (Allen et al. 2013, ApJ 773, 91) demodulates orbits in the time domain with circular-orbit templates, and is sensitive to orbital periods as short as 11 min. Its later CPU→GPU port is mentioned in follow-up papers — [arXiv:1303.0028](https://ar5iv.labs.arxiv.org/html/1303.0028), [arXiv:2007.14849](https://arxiv.org/pdf/2007.14849)

### Inferences
- Choosing a kernel per function:
  - **Isolated and long-period pulsars:** a coherent-harmonic-sum stage (PRESTO `coherentsearch` C code as reference) is now the most sensitive published option and is cheaper than FFA. riptide's FFA is still worth keeping as an independent, well-characterised baseline.
  - **Binaries:** keep FFT accel/jerk with incoherent harmonic summing (PRESTO accelsearch on CPU, PrestoZL or a peasoup-style GPU kernel). Coherent summing has no acceleration dimension yet.
- Published GPU speedups (50–90x) were measured against PRESTO v4/v5. Against v6.3.3 the effective GPU advantage is probably several times smaller, but no one has re-benchmarked. Re-benchmark before investing.
- For a Rust port, the PRESTO v6.3.x accelsearch inner loops are the best-tuned CPU reference: row-wise harmonic summing with cache locality, FMA, and FFTW wisdom. riptide's C++ FFA core (MIT) is the cleanest to port license-wise.

### Gaps
- No single-number FFT+IHS efficiency could be extracted from Morello 2020; it is figure-only.
- No published per-beam peasoup throughput for TRAPUM or MMGPS, and no head-to-head peasoup vs PrestoZL vs AstroAccelerate benchmark on the same data.
- Not verified: the full PrestoZL ApJS paper text (sensitivity tables), or GPU timing details from the Ransom 2026 paper (which GPU, which code).
- Did not read Allen et al. 2013 in full for GPU kernel details.

## Which codes are maintained in 2025–2026, licenses, CPU vs GPU

### Takeaway
Actively maintained in 2026: PRESTO (GPL-2.0, CPU, very active), PrestoZL (GPL-2.0, CUDA), riptide (MIT, CPU), PulsarX (GPL-3.0, CPU), AstroAccelerate (GPL-3.0, CUDA), peasoup (Apache-2.0, CUDA; pushed Aug 2025, with an active vishnubk fork), and gaffa (GPL-3.0, new CUDA FFA). Effectively abandoned: Jintao Luo's presto_on_gpu (2014), Laidler's GPU PRESTO fork (2020), gpu_accelsearch (2017), FFAncy (2017), AMBER (2020, single-pulse/FRB only).

### Cited Findings
| Code | Function | HW | License | Last push | Source |
|---|---|---|---|---|---|
| PRESTO (scottransom) | realfft, rednoise, accelsearch (accel+jerk), coherentsearch, prepfold, ACCEL_sift / sifting.py, coherent_sift.py, stacksearch.py | CPU (OpenMP, FFTW) | GPL-2.0 | 2026-11-18 (v6.3.3 released 2026-11-17) | [repo](https://github.com/scottransom/presto) |
| PrestoZL (zhejianglab) | GPU accelsearch_cu (jerk), prepsubband_cu | CUDA | GPL-2.0 | 2026-08-28 | [repo](https://github.com/zhejianglab/PrestoZL) |
| riptide (v-morello) | FFA, running-median dereddening, rffa pipeline, rseek | CPU (C++ core, Python API) | MIT | 2026-09-21 (v0.2.7, 2026-09-20) | [repo](https://github.com/v-morello/riptide) |
| peasoup (ewanbarr) | dedisp + time-domain resampling accel search + FFT + IHS; now template-bank (circular/elliptical) | CUDA (sm_60–sm_90, CUDA ≤12.6) | Apache-2.0 | 2025-08-05 | [repo](https://github.com/ewanbarr/peasoup) |
| 5D_Peasoup (vishnubk) | template bank | CUDA | GPL-3.0 | retired, merged into peasoup | [repo](https://github.com/vishnubk/5D_Peasoup) |
| AstroAccelerate | dedispersion, FDAS (custom FFT, interbinning, PRESTO-style block-median norm), jerk, harmonic sum, periodicity | CUDA | GPL-3.0 | 2026-09-25 | [repo](https://github.com/AstroAccelerateOrg/astro-accelerate) |
| PulsarX (ypmen) | filtool RFI, dedisperse_all_fil (has `-z zdot`), psrfold_fil folding/optimisation; accel search "under development" | CPU | GPL-3.0 | 2026-11-08 | [repo](https://github.com/ypmen/PulsarX) |
| gaffa (lintian233) | "GPU-accelerated Acceleration-search FFA" package, CUDA/C++ + Python, CUDA 12.8 | CUDA | GPL-3.0 | 2026-08-10 | [repo](https://github.com/lintian233/gaffa) |
| FFAncy (adcameron) | FFA testing package (Cameron et al. 2017) | CPU (C) | MIT | 2017-06-14 | [repo](https://github.com/adcameron/ffancy) |
| ffasearch (vkond) | FFA (Kondratiev) | CPU (C++) | none listed | 2018-01-27 | [repo](https://github.com/vkond/ffasearch) |
| presto_on_gpu (jintaoluo) | GPU accelsearch branch | CUDA | none listed | 2014-05-16 | [repo](https://github.com/jintaoluo/presto_on_gpu) |
| ChrisLaidler/presto | GPU accelsearch fork | CUDA | none listed | 2020-02-08 | [repo](https://github.com/ChrisLaidler/presto) |
| gpu_accelsearch (keithazzopardi) | GPU accel | CUDA | GPL-3.0 | 2017-08-16 | [repo](https://github.com/keithazzopardi/gpu_accelsearch) |
| AMBER (TRASAL) | real-time FRB single-pulse search (not periodicity) | many-core (OpenCL) | Apache-2.0 | 2020-11-17 | [repo](https://github.com/TRASAL/AMBER) |

- PRESTO has no FFA program: a grep of the master tree for `ffa` found nothing. Its periodicity tools are `accelsearch`, `coherentsearch` and `stacksearch.py`. v6.x needs Python ≥3.12, builds with meson, and now pins its regression tests with real data (Ter5A/C/N) — [PRESTO tree](https://github.com/scottransom/presto), [CHANGELOG](https://github.com/scottransom/presto/blob/master/CHANGELOG.md)
- PRESTO v5.1 added `stacksearch.py` (stacking multiple `.fft` files) and `fourier_fold.py` (profiles from complex Fourier amplitudes) — [CHANGELOG](https://github.com/scottransom/presto/blob/master/CHANGELOG.md)
- peasoup's latest features include:
  - template-bank Keplerian searches (ω, τ, φ; plus e, ωp)
  - segmented acceleration search (`--start_sample`, `--fft_size`)
  - `--cdm`, which refines acceleration steps near a coherent DM
  - single-precision harmonic sums
  - PRESTO-compatible `.dat` dumping
  - Docker image `vishnubk/peasoup`

  peasoup does not fold, and recommends PulsarX `psrfold_fil` — [peasoup README](https://github.com/ewanbarr/peasoup)
- 5D_Peasoup's README is marked "retired", with its circular and elliptical template-bank code merged into Peasoup — [5D_Peasoup](https://github.com/vishnubk/5D_Peasoup)
- An active peasoup fork (josephwkania/peasoup, Apache-2.0) was pushed 2026-11-18 — [GitHub search](https://github.com/josephwkania/peasoup)

### Inferences
- License compatibility for a unified codebase:
  - riptide (MIT), peasoup (Apache-2.0) and RustFFT (Apache-2.0/MIT) can be ported or relicensed freely.
  - PRESTO and PrestoZL (GPL-2.0, and PRESTO is not "or later" as far as observed), plus AstroAccelerate, PulsarX and gaffa (GPL-3.0), make derived code GPL. Clean-room reimplementation from the papers is needed if the new stack is permissively licensed.
  - Note that GPL-2.0-only and GPL-3.0 code cannot legally be combined in one binary. Check whether PRESTO's headers say "or later".
- The GPU binary-search choice is peasoup: Apache-2.0, time-domain resampling plus template bank, used in production by TRAPUM/MMGPS. For Fourier-domain jerk it is PrestoZL, which matches PRESTO results exactly.

### Gaps
- Could not list PulsarX's source tree (the API returned 404 for the branch name), so whether PulsarX ships an FFA or accel binary in 2026 is unverified. The README still says accel search is "under development".
- I did not check the GPL-2.0 "or later" wording in PRESTO's headers.
- No Einstein@Home BRP GPU source repository or license was found in this pass.
- Did not inspect gaffa's code or any publication; its README is minimal (CUDA 12.8, Python 3.12, conda).

## Best FFT libraries per platform

### Takeaway
Use FFTW on x86 CPU (PRESTO's whole accelsearch tuning depends on FFTW wisdom), RustFFT/realfft for pure-Rust kernels (SIMD for AVX/SSE4.1/NEON, any size), cuFFT on NVIDIA, and VkFFT for portable GPU (Vulkan/CUDA/HIP/OpenCL/Level Zero/Metal), including Apple GPUs. cuFFT and VkFFT are within about 10% of each other on A100.

### Cited Findings
- FFTW (fftw3) is GPL-2.0, last pushed 2026-10-17 — [FFTW/fftw3](https://github.com/FFTW/fftw3). PRESTO benchmarks are run against FFTW 3.3.10. PRESTO's `makewisdom` defaults to a fast FFTW_MEASURE pass, with `-patient` for FFTW_PATIENT, and v6.3.x asks users to rerun it because wisdom affects speed. It now also plans `coherentsearch`'s FFTs, making them about 1.8x faster — [PRESTO CHANGELOG](https://github.com/scottransom/presto/blob/master/CHANGELOG.md)
- RustFFT (Apache-2.0, pushed 2026-11-21, about 970 stars) is pure Rust and computes any size, including primes, in O(n log n). It picks AVX/FMA, SSE4.1, AArch64 NEON or WASM SIMD automatically (WASM needs the `wasm_simd` feature) — [RustFFT](https://github.com/ejmahler/RustFFT). The companion crate `realfft` (v3.5.1) provides R2C/C2R — [crates.io realfft](https://crates.io/crates/realfft)
- VkFFT (MIT, about 1.8k stars, pushed 2026-04-04) supports Vulkan, CUDA, HIP, OpenCL, Level Zero and Metal backends. It does C2C, R2C/C2R ("up to 2x faster than C2C and 2x less memory") and R2R DCT I–IV, in half, single and double precision — [VkFFT](https://github.com/DTolm/VkFFT)
- Independent TurboFFT study on A100: VkFFT averages 10% overhead against cuFFT (range +3% to −51% by case), and 11% in double precision — [arXiv:2405.02520](https://arxiv.org/pdf/2405.02520). The VkFFT SC22 poster claims better radix (2–13), Bluestein and Rader coverage than cuFFT/rocFFT, and notes rocFFT lacks Rader — [SC22 poster](https://sc22.supercomputing.org/proceedings/tech_poster/tech_poster_pages/rpost143.html)
- cuFFT is efficient for sizes whose prime factors are 2, 3, 5 or 7, but performance varies by GPU. peasoup advises benchmarking FFT sizes and setting `--fft_size` explicitly — [peasoup README](https://github.com/ewanbarr/peasoup)
- AstroAccelerate's FDAS has an option to use its own custom FFT (`fdas_custom_fft`) instead of the library FFT — [AstroAccelerate README](https://github.com/AstroAccelerateOrg/astro-accelerate)

### Inferences
- Pulsar FFTs are large 1-D real transforms (2^23–2^27 points; peasoup's example uses fft_size 67,108,864). Accel/jerk search then needs many small, batched complex FFTs or IFFTs for the f-fdot correlation, typically lengths of a few thousand. Choose libraries per call pattern:
  - **Big R2C:** FFTW on CPU, cuFFT on NVIDIA.
  - **Batched small C2C:** cuFFT batched, VkFFT, or fused custom kernels as in AstroAccelerate and Dimoudi.
- For a Rust CPU kernel, RustFFT avoids FFTW's GPL. Its speed relative to FFTW for 2^24+ real transforms was not verified here.
- For Apple Silicon, VkFFT's Metal backend is the only cross-vendor GPU option found. CPU options are Accelerate/vDSP or RustFFT with NEON.

### Gaps
- No primary-source benchmark found comparing RustFFT and FFTW at pulsar-scale lengths, or for Apple vDSP/Accelerate FFT. vDSP relevance is therefore inferred, not cited.
- cuFFT licensing (proprietary, redistributable runtime) and cuFFTDx device-callable FFTs were not researched.

## Any Rust implementations of FFA / accel search?

### Takeaway
None found. GitHub and crates.io searches turned up no Rust pulsar periodicity-search code (FFA, accelsearch, harmonic summing or dedispersion). "pulsar" in Rust is dominated by Apache Pulsar messaging clients. A Rust implementation would be new work. The radio-astronomy Rust ecosystem exists but covers interferometry and reduction, not pulsar search.

### Cited Findings
- `gh search repos "pulsar" --language rust` returned only Apache Pulsar clients and a game engine, and a `dedispersion` Rust search returned nothing. Rust radio-astronomy hits are rubbl-rxpackage (reduction utilities), casacure (casacore replacement) and fringe (simulator) — [GitHub search](https://github.com/search?q=pulsar+language%3ARust&type=repositories), [rubbl-rxpackage](https://github.com/pkgw/rubbl-rxpackage), [casacure](https://github.com/tmolteno/casacure)
- A crates.io query for "pulsar fold" returned only Apache Pulsar client crates (the magnetar-* family) and a VLBI fringe-search crate (frinZ) — [crates.io search](https://crates.io/search?q=pulsar%20fold)
- Non-Python/C FFA reimplementations do exist, e.g. an Ada FFA (MIT, 2026) — [Ada-Fast-Folding-Algorithm](https://github.com/RobertBoettcherSF/Ada-Fast-Folding-Algorithm)

### Inferences
- The best porting references for Rust CPU kernels:
  - **FFA:** riptide's C++ core (MIT; `periodogram`, running median, `rollback`, boxcar S/N with variance-corrected downsampling — note the v0.2.6 fix for downsampled variance near-integer factors) — [riptide CHANGELOG](https://github.com/v-morello/riptide/blob/master/CHANGELOG.md)
  - **Coherent harmonic sum:** the PRESTO `coherentsearch` C and CoherentSearch.jl. The Julia code's license should be checked; PRESTO's port is GPL.
  - **Accel/jerk with IHS:** PRESTO v6.3.3 accelsearch (GPL, algorithm reference only for a permissive stack).

### Gaps
- Search coverage was limited to GitHub repository search and crates.io keyword search. Private or unlisted Rust code, and code without "pulsar" in its description, could exist.

## Candidate sifting (ACCEL_sift.py, coherent_sift.py, peasoup distillers)

### Takeaway
PRESTO's `sifting.py` / `ACCEL_sift.py` is still maintained, with fixes in 2026, and PRESTO now also ships `coherent_sift.py`, which ranks candidates by pulsar-likeness and flags harmonics, known pulsars and RFI. peasoup does not fold; TRAPUM-style pipelines feed peasoup XML candidates to PulsarX for folding.

### Cited Findings
- `examplescripts/ACCEL_sift.py` now gets DMs from `sifting.DM_re`; `sifting`'s "loglin" axis was fixed and tests were added (`tests/test_sifting.py`) — [PRESTO CHANGELOG](https://github.com/scottransom/presto/blob/master/CHANGELOG.md), [PRESTO tree](https://github.com/scottransom/presto)
- `coherent_sift.py` sifts `.cohout` files across DMs and observations. It groups detections, flags harmonics, known pulsars, deep-pass signals and edge-peaked or broadband (RFI-like) candidates, ranks by pulsar-likeness, and writes HTML summaries — [PRESTO CHANGELOG v6.2.0](https://github.com/scottransom/presto/blob/master/CHANGELOG.md)
- PRESTO v5.0 added `compare_periods.py`, which matches candidate periods and harmonics against parfiles — [CHANGELOG](https://github.com/scottransom/presto/blob/master/CHANGELOG.md)
- riptide detects candidates with a dynamic threshold: frequency segments with median and IQR-based sigma, fit by a polynomial in log(ν) — [ar5iv 2004.03701](https://ar5iv.labs.arxiv.org/html/2004.03701)
- PrestoZL's group published Pulsar-VLM (ApJS), a multimodal model that classifies PRESTO/PrestoZL diagnostic plots with 94.2–100% recall across five telescopes — [PrestoZL README](https://github.com/zhejianglab/PrestoZL)

### Inferences
- A unified stack's sifting layer should be pure Python orchestration: DM clustering plus harmonic matching, as in ACCEL_sift, extended with coherent_sift's RFI and known-pulsar flags. It should consume a common candidate schema emitted by all search backends (accel, coherent, FFA, template bank).

### Gaps
- Did not retrieve peasoup's candidate "distiller" code (the harmonic, DM and acceleration distillers in its C++ source) or its exact merging tolerances. riptide's internal harmonic-flagging rules were not read.
