# Dedispersion and Single-Pulse / FRB Transient Search Implementations (state as of Nov 2026)

Scope: incoherent dedispersion (brute force, tree, subband, FDMT, Fourier-domain), DM planning, single-pulse detection (boxcar matched filter + candidate clustering), coherent dedispersion; for picking best-of-breed kernels for a unified Rust CPU / C++-CUDA GPU / Python orchestration stack. Repo metadata (license, last push, stars) was pulled from the GitHub API on 2026-11-24.

## Published benchmarks: which is fastest on CPU and on GPU?

### Takeaway
The newest and best-documented cross-engine numbers come from `dmt` (pravirkr, MIT, C++20 + CUDA/HIP + Python): FDMT at 55x real time on 8 EPYC cores and 776x real time on one L40S for a MeerKAT L-band setup (4096 ch, ~4057 DMs), ~19x (CPU) / ~6x (GPU) faster than its own brute-force engine. For brute force on GPU, AstroAccelerate (Novotny et al. 2023) claims ~4x over other available solutions; TransientX is the reference CPU subband single-pulse searcher (real time on one core for 4096 ch x 4096 DMs at 153 us).

### Cited Findings
**GPU brute force / classic codes**
- Barsdell et al. 2012 (dedisp library, basis of heimdall): direct GPU dedispersion on GTX 480 ~9x faster than an optimised quad-core i7 930; Tesla C1060 ~5x; a single GPU computes the transform in < 1/3 of observation time for 1 min of 2-bit data, 1024 ch, 64 us, 1196 DMs (0-1000). CPU direct scaled only ~2.5x on 4 cores (memory-bandwidth bound). Their direct code ~2.3x faster than Magro et al. 2011 on the same C1060. — [Barsdell et al. 2012, ar5iv](https://ar5iv.labs.arxiv.org/html/1201.5380)
- Barsdell 2012 tree/subband speedups over direct: piecewise-linear tree up to ~3x (N'nu=64); subband up to ~20x (N'nu=256) in discussion, while abstract/conclusions quote 3-7x for both (internally inconsistent figures); frequency-padded tree ~1.5x. — [Barsdell et al. 2012](https://ar5iv.labs.arxiv.org/html/1201.5380)
- AstroAccelerate brute-force dedispersion (Novotny, Adamek, Clark, Giles, Armour; ApJS 269, 2023): ~4x faster than other available solutions for 8-bit, 8192-channel input; >=20x faster than real time for 11 existing telescope configurations; runs on NVIDIA GPUs and Intel/AMD CPUs. — [arXiv:2311.05341](https://arxiv.org/abs/2311.05341)
- Heimdall++ (Xia et al., Nov 2025): re-engineered heimdall with fine-grained GPU parallelism, better memory management, CPU/GPU stage separation; up to 2.66x faster single-file and 2.05x multi-file batch vs original heimdall on RTX 3080 Ti, results consistent with original. No repo/license stated. — [arXiv:2512.00398](https://arxiv.org/abs/2512.00398)
- Heimdall validation on Northern Cross (Cesare et al., Sep 2026): one beam per GPU; real-time PSRDADA mode "stays in the real-time flow for most runs"; PSRDADA vs SIGPROC modes agree within 3% in TOA/DM, 10% in S/N. No numeric real-time factor. — [arXiv:2609.20635](https://arxiv.org/abs/2609.20635)
- ASTROFLOW (Nov 2025; ApJS DOI 10.3847/1538-4365/ae4a26): C++/CUDA core with Python interface, shared-memory tiling, OpenMP; ~10x faster than real time on a 150 s, 2048-channel observation on consumer GPUs; YOLO-based detector on dynamic-spectrum tiles. — [arXiv:2511.02328](https://arxiv.org/abs/2511.02328); repo BSD-3-Clause, pip package `pulseflow`, last push 2026-06 — [lintian233/astroflow](https://github.com/lintian233/astroflow)
- ASTROFLOW paper reports heimdall and TransientX recover most astrophysical bursts on FAST data but produce substantially more false detections (than ASTROFLOW's ML detector). — [arXiv:2511.02328 (search snippet)](https://arxiv.org/pdf/2511.02328)

**FDMT / tree**
- `dmt` v0.7.0 benchmark: AMD EPYC 9654 (8 threads, one per CCD) FDMT 55x real time; NVIDIA L40S (device-resident) 776x real time; ~19x faster than brute-force DDMT on CPU and ~6x on GPU, `meerkat_lband` preset (4096 ch, 856-1712 MHz, 612.48 us, 20.1 s streaming blocks, up to 4057 DMs). — [pravirkr/dmt README](https://github.com/pravirkr/dmt)
- Zackay & Ofek FDMT: exact incoherent dedispersion with complexity N_f N_t log2(N_f); authors claim a single-CPU implementation is competitive with brute-force GPU dedispersion (benchmark N_f=2^10, N_t=5x2^16, N_DM=2^10). — [arXiv:1411.5373](https://arxiv.org/pdf/1411.5373)
- FREDDA (ASKAP/CRAFT GPU FDMT pipeline): blocks of 512 integrations per beam, FDMT computing 4096 DM trials at 0.92 pc/cc resolution, then boxcar convolution per DM trial. — [Bannister et al. 2017, ApJL](https://iopscience.iop.org/article/10.3847/2041-8213/aa71ff)
- Bifrost stream framework ships a GPU FDMT block (`bf.blocks.fdmt`); BSD-3-Clause, actively pushed (2026-11). — [lwa-project/bifrost](https://github.com/lwa-project/bifrost)
- Julia FastDMTransform.jl: CPU FDMT covers all DMs in a 3 s chunk in ~50 ms (~10x faster than reference Python); its CUDA comparison (55-72 ms) is brute force, not FDMT on GPU. — [FastDMTransform.jl](https://github.com/kiranshila/FastDMTransform.jl)

**CPU single-pulse**
- TransientX (Men & Barr 2024): CPU subband dedispersion with AVX2 + OpenMP; real time on a single CPU core for 4096 DMs, 4096 channels, 153 us (Intel i7-10750H, simulated 600 s MMGPS-L-like data); ~an order of magnitude faster than BEAR; no numeric PRESTO/heimdall comparison. — [arXiv:2401.13834](https://arxiv.org/html/2401.13834v1)
- `sbirat`: Python (nanobind) binding of the Novotny et al. 2023 fast CPU brute-force algorithm (original C at jan2nov/dedispersion); MIT. — [astrogewgaw/sbirat](https://github.com/astrogewgaw/sbirat)

**Fourier-domain**
- Fourier-domain dedispersion (FDD, Bassa/ASTRON line of work) reduces dedispersion time vs GPU brute force but at the cost of longer preprocessing. — [arXiv:2110.03482](https://arxiv.org/pdf/2110.03482); ASTRON fork of dedisp used for the FDD project — [git.astron.nl/RD/dedisp](https://git.astron.nl/RD/dedisp)
- Direct dedispersion on many-core accelerators is inherently memory-bandwidth-bound (Sclocco et al. AMBER line). — [arXiv:1601.05052](https://arxiv.org/html/1601.05052)

### Inferences
- GPU fastest: FDMT-on-GPU (`dmt`, ~776x RT on L40S) beats brute force by ~6x on the same hardware; among brute-force codes AstroAccelerate is the published leader. Heimdall's dedisp is older and slower than both (Heimdall++ shows 2-2.7x headroom just from engineering).
- CPU fastest: FDMT (`dmt` CPU, 55x RT on 8 EPYC threads) for dense DM grids; TransientX subband is close in practice due to cache friendliness and is a complete pipeline. Brute force on CPU is memory-bound (~2.5x scaling on 4 cores in 2012).
- For a unified stack, `dmt` (MIT, C++/CUDA/HIP, Python bindings, streaming, packed 1-16 bit input, per-trial calibrated noise) is the strongest candidate to wrap or port for the core transform; its algorithms are the natural reference for a Rust CPU port.
- Benchmarks across papers are not apples-to-apples (different N_chan, N_DM, bit depth, tsamp, GPU generation); a common harness is needed before final selection.

### Gaps
- No peer-reviewed head-to-head of dedisp vs AstroAccelerate vs FDMT-GPU vs TransientX on identical data/hardware was found.
- AMBER/TRASAL absolute throughput numbers not retrieved. CHIME Bonsai throughput numbers not retrieved (repo not found under kmsmith137/bonsai or CHIMEFRB/bonsai via API).
- DM-trials/s numbers for heimdall in 2025-2026 hardware not found (only relative Heimdall++ speedups).
- Absolute TransientX matched-filter/clustering timings exist only in figures of the paper (not extracted).

## Which are actively maintained in 2025-2026? Licenses?

### Takeaway
Actively maintained (pushes in H2 2026): PRESTO (GPL-2.0), TransientX/PulsarX (GPL-3.0), AstroAccelerate (GPL-3.0), sigpyproc3 (MIT), dmt (MIT), bifrost (BSD-3), cdmt (no license file), pirate (no license file), astroflow (BSD-3). Permissive licenses (MIT/BSD) are only dmt, sigpyproc3, bifrost, astroflow, PyTorchDedispersion, sbirat, chime-sps/FDMT; the main CPU/GPU search engines (PRESTO, TransientX, AstroAccelerate) are GPL, which constrains code reuse in a permissively licensed unified codebase.

### Cited Findings
(License / last push / stars from GitHub API, 2026-11-24)
- pravirkr/dmt — MIT, pushed 2026-11-21, 4 stars; C++20, CUDA 12.2+ and HIP, Python 3.12+, CI + codecov, ReadTheDocs; engines FDMT, FDMT-FFT, DDMT (brute), SDMT, DDMT-FFT, CohFDMT (coherent search of baseband). — [github.com/pravirkr/dmt](https://github.com/pravirkr/dmt)
- scottransom/presto — GPL-2.0, pushed 2026-11-18, 290 stars. — [github.com/scottransom/presto](https://github.com/scottransom/presto)
- ypmen/TransientX — GPL-3.0, pushed 2026-07-24; ypmen/PulsarX — GPL-3.0, pushed 2026-11-08. — [TransientX](https://github.com/ypmen/TransientX), [PulsarX](https://github.com/ypmen/PulsarX)
- AstroAccelerateOrg/astro-accelerate — GPL-3.0, pushed 2026-09-25; used in MeerTRAP. — [repo](https://github.com/AstroAccelerateOrg/astro-accelerate); [arXiv:2311.05341](https://arxiv.org/abs/2311.05341)
- FRBs/sigpyproc3 — MIT, pushed 2026-10-17. — [repo](https://github.com/FRBs/sigpyproc3)
- lwa-project/bifrost (GPU FDMT block) — BSD-3-Clause, pushed 2026-11-18. — [repo](https://github.com/lwa-project/bifrost)
- lintian233/astroflow — BSD-3-Clause, pushed 2026-06-26, CUDA 12.6+, PyPI `pulseflow`. — [repo](https://github.com/lintian233/astroflow)
- kmsmith137/pirate ("Perimeter Institute RAdio Transient Engine") — experimental GPU fast-transient search for CHORD (successor to CHIME Bonsai line by K. Smith); real-time server code + GPU dedispersion, currently no RFI flagging, Ubuntu + CUDA only; no license file; pushed 2026-10-28. — [repo](https://github.com/kmsmith137/pirate)
- cbassa/cdmt (Coherent Dispersion Measure Trials) — no license file in API, pushed 2026-10-03. — [repo](https://github.com/cbassa/cdmt); paper — [arXiv:1607.00909](https://arxiv.org/pdf/1607.00909)
- ewanbarr/dedisp — multi-GPU fork of Barsdell dedisp, "works with CUDA > 12", no license file, pushed 2025-04. ajameson/dedisp — no license file, pushed 2024-06. — [ewanbarr/dedisp](https://github.com/ewanbarr/dedisp), [ajameson/dedisp](https://github.com/ajameson/dedisp)
- nkosogor/PyTorchDedispersion — BSD-3-Clause, pushed 2025-08/11; incoherent dedispersion + boxcar + candidate detection in PyTorch; used for OVRO-LWA Time Machine searches. — [repo](https://github.com/nkosogor/PyTorchDedispersion); [arXiv:2512.18543](https://arxiv.org/pdf/2512.18543)
- DSA-110 "cerberus" search uses JAX JIT dedispersion + boxcar adapted from PyTorchDedispersion, precomputing per-channel shift indices on GPU. — [arXiv:2510.18136](https://arxiv.org/pdf/2510.18136)
- chime-sps/FDMT — Python FDMT, MIT, pushed 2025-08. — [repo](https://github.com/chime-sps/FDMT)
- TRASAL/AMBER — Apache-2.0, last push 2020-11 (stale). — [repo](https://github.com/TRASAL/AMBER)
- jintaoluo/presto_on_gpu — no license, last push 2014 (dead). — [repo](https://github.com/jintaoluo/presto_on_gpu)
- activexray/FastDMTransform.jl — MIT, 2022-2024. Julia family (astrogewgaw/Dedispersion.jl, FastDedispersion.jl, DirectDedispersion.jl, activexray/Dedisp.jl) last updated 2024-2025. — [GitHub search results](https://github.com/astrogewgaw/FastDedispersion.jl)
- CHIME/Pulsar semi-coherent dedispersion search pipeline (kmjc/chime_psr_survey_public), pushed 2026-01. — [repo](https://github.com/kmjc/chime_psr_survey_public)

### Inferences
- For a permissively licensed unified codebase, `dmt` (MIT) is the only modern, maintained, multi-backend dedispersion core; GPL codes (PRESTO, TransientX, AstroAccelerate) can be used as references/benchmarks or invoked as external tools but copying code imposes GPL on the result.
- dedisp/heimdall lineage is in maintenance-fork mode (personal forks, no explicit license files), so it is a poor foundation for new code.
- AMBER is effectively unmaintained since 2020.

### Gaps
- `your` (thepetabyte/your) and heimdall canonical repo locations could not be resolved via GitHub API (404 on guessed names); their 2026 maintenance and license status is unverified here (heimdall historically on SourceForge).
- DSPSR license/maintenance not checked in this pass.
- License for dedisp is not stated in the repos checked.

## Algorithmic trade-offs (accuracy / S/N loss of tree / FDMT vs brute force; subband smearing)

### Takeaway
Brute force is exact but O(N_DM N_t N_f) and memory-bound; FDMT is exact (with box-smearing handling) at O(N_t N_f log N_f) but yields a fixed delay grid set by tsamp; subband dedispersion is O(N_DM N_t sqrt(N_f)) and cache-friendly but adds smearing that grows quickly with subband size; the piecewise-linear tree is faster but loses S/N from curvature approximation.

### Cited Findings
- Barsdell 2012: piecewise-linear tree adds <10% extra smearing at DM 1000 for subband size <=64 (64 judged optimal); subband method adds ~40% extra smearing at N'nu=16 and ~300% at N'nu=32, so sizes >32 unlikely acceptable; frequency-padded tree adds no extra smearing but needs more channels. Conclusion: GPU direct fast enough for real-time; tree ~3x faster at modest signal loss but complex; subband fast and simple but smearing substantial; no single best algorithm. — [Barsdell et al. 2012](https://ar5iv.labs.arxiv.org/html/1201.5380)
- TransientX: subband complexity O(N_DM N_t sqrt(N_f)) vs O(N_DM N_t log2 N_f) for tree/FDMT; speed "comparable to FDMT" because more cache friendly; supports DM plans that downsample before dedispersion per DM range to address intra-channel smearing. — [arXiv:2401.13834](https://arxiv.org/html/2401.13834v1)
- TransientX S/N losses: spectral index 2 with eta=2 gives ~6% loss; boxcar vs scattered Gaussian ~9% even at scattering 10x width. — [arXiv:2401.13834](https://arxiv.org/html/2401.13834v1)
- FDMT claims exact incoherent transform at N_f N_t log2 N_f. — [Zackay & Ofek, arXiv:1411.5373](https://arxiv.org/pdf/1411.5373)
- `dmt` exposes `use_box_smearing` in FDMT and computes exact per-trial noise sigma so S/N needs no calibration run; all engines share plans/DM grids so they can be swapped and compared on same data. — [pravirkr/dmt](https://github.com/pravirkr/dmt)
- Fourier-domain dedispersion: arbitrary DM grids with less dedispersion time but larger preprocessing cost. — [arXiv:2110.03482](https://arxiv.org/pdf/2110.03482)
- pFDMT (pruned, cache-friendly FDMT) targets candidate folding with broadly distributed DMs rather than single-pulse search. — [arXiv:2309.02544](https://arxiv.org/pdf/2309.02544)
- BASSET explores sub-band pulse search strategies for faint narrow-band FRBs (GPU dedispersion + OpenMP). — [arXiv:2501.05875](https://arxiv.org/pdf/2501.05875)
- STRIDE (Di Pietrantonio et al., Mar 2026): streaming dedispersion for imaging-based (per-pixel) searches, cutting minimum memory from 684.5 GB to 14.4 GB in the test case; no code link. — [arXiv:2603.02931](https://arxiv.org/abs/2603.02931)

### Inferences
- Recommended design: brute-force (exact, arbitrary DM grid) as the reference/validation engine; FDMT as the production engine for dense DM grids; subband (PRESTO DDplan-style or TransientX-style with per-range downsampling) as CPU fallback for very high channel counts. Expose all behind a single plan object as `dmt` does.
- Narrow-band FRBs argue for sub-band searches (BASSET) as an optional mode.

### Gaps
- Quantitative S/N loss of FDMT vs brute force across DM (due to delay-grid quantization at the tsamp lattice) was not found in a peer-reviewed comparison in this pass.
- PRESTO DDplan.py algorithm details (smearing budget combining tsamp, channel, DM-step, subband terms) were not re-fetched; heimdall's DM-list generation (pulse-width tolerance, typical 1.25) not verified from source here.

## Single-pulse detection and candidate clustering

### Takeaway
All modern codes use boxcar matched filtering over a geometric set of widths on each dedispersed series, then cluster in (DM, time, width): TransientX uses DBSCAN with a kd-tree; heimdall uses its own GPU candidate merging; newer pipelines (ASTROFLOW, DSA-110) add ML detectors to cut false positives.

### Cited Findings
- TransientX: boxcar widths in a geometric sequence with tunable ratio, O(N_t) per series per width; DBSCAN in C++ with kd-tree neighbor search (O(log N) per lookup) on the DM-time plane after collapsing width by max S/N, DM converted to dispersion delay; benchmark radii 4, 8, 16 samples on 1 s blocks; clustering time linear in number of points. Companion plotting library PlotX. — [arXiv:2401.13834](https://arxiv.org/html/2401.13834v1)
- FREDDA: boxcar convolution per DM trial after FDMT. — [Bannister et al. 2017](https://iopscience.iop.org/article/10.3847/2041-8213/aa71ff)
- PyTorchDedispersion: incoherent dedispersion + boxcar range + candidate detection in one GPU pipeline. — [arXiv:2512.18543](https://arxiv.org/pdf/2512.18543)
- Heimdall transfers only candidate metadata (arrival time, DM, S/N) back to host in the Northern Cross pipeline. — [search snippet, Northern Cross](https://arxiv.org/pdf/2609.20635)
- GSP (FAST/CRAFTS) GPU single-pulse pipeline with database uses PRESTO-style subband dedispersion on GPU. — [arXiv:2110.12749](https://arxiv.org/pdf/2110.12749)
- ASTROFLOW replaces threshold+cluster with a YOLO detector on dynamic-spectrum tiles. — [arXiv:2511.02328](https://arxiv.org/abs/2511.02328)

### Inferences
- A Rust CPU boxcar via prefix sums (O(N_t) per width) plus DBSCAN/FoF on (DM-delay, time) mirrors TransientX's proven design and avoids GPL copying if reimplemented from the paper.
- Calibrated per-trial noise (as in `dmt`) removes the need for running-median normalization hacks in S/N estimation.

### Gaps
- PRESTO single_pulse_search.py details (FFT-convolution boxcars, downfacts list, sigma threshold), AstroAccelerate SPDT, and heimdall's merge/"giant" clustering rules were not verified from primary sources in this pass. "presto-sp/spcand" not located.

## Rust, Apple Silicon/Metal, and portable GPU (HIP/SYCL/Vulkan) implementations

### Takeaway
No Rust dedispersion/FDMT crate exists on crates.io or GitHub (only a Rust SIGPROC filterbank reader). HIP is covered by `dmt`; SYCL by a coherent-dedispersion backend (fxzjshm); Metal by one tiny hobby repo (dedispulse). Vulkan compute: none found.

### Cited Findings
- crates.io searches for "dedispersion" and "fdmt" returned no relevant crates (only an unrelated VLBI fringe-search crate); GitHub Rust searches for pulsar/radio transient return nothing relevant. Only relevant Rust code: `sigproc_filterbank` (pure-Rust SIGPROC filterbank read/write, 2024) and GReX `dada2fil` (PSRDADA to filterbank, 2025). — [activexray/sigproc_filterbank](https://github.com/activexray/sigproc_filterbank), [GReX-Telescope/dada2fil](https://github.com/GReX-Telescope/dada2fil)
- `dmt`: CPU SIMD AVX2/AVX-512/NEON (so ARM / Apple Silicon CPU supported), GPU CUDA and HIP with identical results to float rounding; no Metal. — [pravirkr/dmt](https://github.com/pravirkr/dmt)
- eigenpaul/dedispulse: Apple Silicon Metal incoherent dedispersion + folding for SIGPROC .fil (validated on CAMRAS B0329+54), single-DM/sweep tool, no license, 0 stars, May 2026. — [repo](https://github.com/eigenpaul/dedispulse)
- fxzjshm/simple-radio-telescope-backend: SYCL 2020 (AdaptiveCpp, oneAPI DPC++) real-time coherent dedispersion + RFI + single-pulse detection from baseband; tested on CUDA and ROCm/HIP; OpenCL/SPIR-V and Vulkan (clspv+clvk) paths untested. — [repo](https://github.com/fxzjshm/simple-radio-telescope-backend); [FPS12 slides, Zhou](https://psr.pku.edu.cn/conference/fps/FPS12/slide/zhou_yukai.pdf)
- Earlier portable OpenCL dedispersion (Sclocco et al. 2016) tested on AMD/NVIDIA GPUs and Xeon Phi. — [arXiv:1601.05052](https://arxiv.org/html/1601.05052); sievers/opencl_dm_transform (2015) — [repo](https://github.com/sievers/opencl_dm_transform)
- AstroAccelerate brute force runs on Intel/AMD CPUs as well as NVIDIA GPUs. — [arXiv:2311.05341](https://arxiv.org/abs/2311.05341)

### Inferences
- A Rust FDMT/brute-force/boxcar CPU kernel set would be novel; port from `dmt` (MIT) or the Zackay & Ofek paper, with `std::simd`/NEON for Apple Silicon.
- For Apple GPUs, Metal (via `metal-rs` / wgpu) would be greenfield; wgpu/WGSL compute could give Metal + Vulkan + DX12 portability from Rust.
- HIP support exists today via `dmt`; SYCL only for coherent dedispersion.

### Gaps
- No benchmarks of any Metal, Vulkan, or SYCL incoherent dedispersion were found.

## Coherent dedispersion (DSPSR, cdmt, GPU coherent codes)

### Takeaway
cdmt (Bassa et al.) remains the reference for GPU coherent DM trials for searches and is still being pushed in 2026; `dmt` adds CohFDMT for coherent search of baseband; SYCL backend by fxzjshm offers vendor-neutral real-time coherent dedispersion.

### Cited Findings
- cdmt ("Coherent Dispersion Measure Trials"), GPU coherent dedispersion over many DM trials enabling pulsar/fast-transient searches (LOFAR). — [arXiv:1607.00909](https://arxiv.org/pdf/1607.00909); [cbassa/cdmt](https://github.com/cbassa/cdmt) (pushed 2026-10-03, no license file)
- `dmt` CohFDMT: coherent search of baseband voltages within the same library. — [pravirkr/dmt](https://github.com/pravirkr/dmt)
- GMRT real-time coherent dedispersion pipeline (GPU). — [arXiv:1509.00186](https://arxiv.org/pdf/1509.00186)
- CoherentDedispersion.jl (Julia, david-macmahon, 2026-07); ypmen/BasebandX (coherent dedispersion in development, folding available, 2025-12). — [CoherentDedispersion.jl](https://github.com/david-macmahon/CoherentDedispersion.jl), [BasebandX](https://github.com/ypmen/BasebandX)

### Inferences
- For the unified stack, coherent dedispersion is lower priority; wrap cdmt/DSPSR externally, and consider `dmt` CohFDMT as the MIT-licensed in-library option.

### Gaps
- DSPSR GPU performance numbers and current maintenance/license not fetched. David McKenna's cdmt fork (David-McKenna/cdmt) last pushed 2021.
