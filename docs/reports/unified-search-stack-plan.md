# Build a permissive core, wrap the GPL giants

The fastest route to a unified, fast pulsar and transient codebase is not to port PRESTO, DSPSR and PSRCHIVE into one program. Build a **small, MIT-licensed core** instead: one `definitions/` folder holding headers, DM plans and candidate schemas; Rust CPU kernels; C++/CUDA GPU kernels; and a Python orchestrator that picks a backend at runtime. Then **wrap or call the established codes as external tools and test oracles**, rather than copying their source. Licensing drives this more than taste. The engines that dominate each stage (PRESTO, PrestoZL, TransientX, PulsarX, AstroAccelerate, `your`, FETCH) are GPL. The permissive codes that are both best-in-class and maintained are few: `dmt` for dedispersion, riptide for the FFA, peasoup for GPU binary search, sigpyproc3 for I/O and IQRM, and clfd for cleaning folded data. Those few are what you build on. For your situation, Mac development with no CUDA and a real-time transient pipeline for the Diamond Array Telescope plus a two-element interferometer, the order matters. Do the single-pulse chain first (I/O, RFI, dedispersion, boxcar search, clustering, ML vetting) with the CPU backend as the reference implementation. Add CUDA on a Linux GPU host second, and periodicity search third. The repository change you asked for is a precondition: rename `psr-searchmode-to-filterbank` to `filterbank` and move the converter into `filterbank/converter/`. The one decision only you can make is whether `filterbank` stays a git submodule. I recommend converting it into a regular directory, keeping the converter's history.

## The licence line, not the language line, sets the architecture

Most of the best search code is GPL. **PRESTO is GPL-2.0** ([scottransom/presto](https://github.com/scottransom/presto)). TransientX and PulsarX are GPL-3.0 ([TransientX](https://github.com/ypmen/TransientX), [PulsarX](https://github.com/ypmen/PulsarX)), and so are AstroAccelerate ([repo](https://github.com/AstroAccelerateOrg/astro-accelerate)) and the GPU accelsearch fork PrestoZL, which is GPL-2.0 ([PrestoZL](https://github.com/zhejianglab/PrestoZL)). On the Python side, `your` and FETCH are both GPL-3.0 ([your](https://github.com/thepetabyteproject/your), [FETCH](https://github.com/devanshkv/fetch)). PSRCHIVE and DSPSR are listed under the Academic Free License ([psrchive](https://sourceforge.net/projects/psrchive/), [dspsr](https://sourceforge.net/projects/dspsr/)). The AFL is generally regarded as GPL-incompatible, so AFL and GPL code should not be linked into one binary.

A line-by-line translation of PRESTO's C into Rust is still a modified version of PRESTO, so the result must be GPL-2.0. The FSF FAQ addresses translation to another language directly ([GNU GPL FAQ](https://www.gnu.org/licenses/gpl-faq.html)), although its exact wording could not be re-fetched for this report. There is a further trap: GPL-2.0-only and GPL-3.0 code cannot be combined in one binary. A stack that "ports the best of everything" could therefore end up legally unbuildable, even before any question of permissiveness.

The permissive codes cover more ground than you might expect, and they are the modern ones. **`dmt` is MIT**, actively pushed in November 2026, written in C++20 with CUDA, HIP and Python bindings, and it ships FDMT, brute-force, subband and coherent engines ([pravirkr/dmt](https://github.com/pravirkr/dmt)). riptide (FFA) and sigpyproc3 (I/O and IQRM) are MIT ([riptide](https://github.com/v-morello/riptide), [sigpyproc3](https://github.com/FRBs/sigpyproc3)). **peasoup is Apache-2.0** and is in production on MeerKAT/TRAPUM ([peasoup](https://github.com/ewanbarr/peasoup)). ASTROFLOW and bifrost are BSD-3 ([astroflow](https://github.com/lintian233/astroflow), [bifrost](https://github.com/lwa-project/bifrost)). clfd is MIT ([clfd](https://github.com/v-morello/clfd)). RustFFT is Apache-2.0/MIT and VkFFT is MIT ([RustFFT](https://github.com/ejmahler/RustFFT), [VkFFT](https://github.com/DTolm/VkFFT)).

The recommended strategy has three tiers:

1. **Core.** License the core MIT. MIT rather than Apache-2.0, because Apache-2.0 is incompatible with GPLv2-only code, and an MIT core stays compatible with a GPL tool someone later chooses to link.
2. **Port.** Freely port or vendor algorithms from the permissive codes listed above.
3. **Call.** Treat every GPL or AFL tool as a **subprocess adapter or test oracle**: run `prepsubband`, `accelsearch`, `psrfold_fil`, `pat` or `digifil`, and read their output files. Reading another program's file format does not make your code a derivative work.

If you later want a true port of a GPL kernel (for example PRESTO's v6.3.3 accelsearch inner loops), put it in a separately distributed, GPL-licensed optional package, so the core wheel stays clean. Keep a `PROVENANCE.md` in each Rust crate and CUDA module stating "clean-room from paper X" or "ported from MIT project Y at commit Z".

Two follow-ups are needed before you commit to this:

- **Check licences that could not be confirmed.** heimdall's licence, the CoherentSearch.jl licence, and whether PRESTO's headers say "or later" were not confirmed in the research notes. Check them before relying on any of these as more than an oracle.
- **Choose a licence for the converter.** Your own converter code needs an explicit licence too. MIT is the natural fit.

## One definitions folder, three languages, zero-copy boundaries

"Definitions in one place, everything else independent" works if each module is a pure function over **arrays plus Arrow tables**. Every module reads its types from a single schema source, and every language boundary uses a standard zero-copy protocol.

**Bulk data.** Dynamic spectra and time series cross language boundaries as plain N-d arrays through DLPack or the buffer protocol. On the Python–C++ side, nanobind's `nb::ndarray` exchanges arrays zero-copy with NumPy, CuPy, PyTorch, JAX and MLX, and can constrain dtype, device and layout ([nanobind ndarray](https://nanobind.readthedocs.io/en/latest/ndarray.html)). On the Python–Rust side, PyO3 with rust-numpy fills the same role.

**Records.** Candidate lists, DM-trial metadata and sifted lists are **Arrow tables**. The Arrow PyCapsule Interface lets pyarrow, polars and Rust `arrow` exchange tables without copying ([Arrow PyCapsule Interface](https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html), [pyo3-arrow](https://docs.rs/crate/pyo3-arrow/latest)).

**Small structs and codegen.** Filterbank header, DM plan and fold parameters are defined once and code-generated into Rust (serde structs), Python (dataclasses or msgspec), C++ (POD headers) and Arrow schema JSON. The best precedent is Rerun: it defines its types in one IDL crate, generates the Rust, Python and C++ SDKs from it, and uses Arrow for all data ([rerun ARCHITECTURE.md](https://raw.githubusercontent.com/rerun-io/rerun/main/ARCHITECTURE.md)). CI should fail if regenerating the code produces a diff.

**Shared numerics.** One design choice removes the most common cross-backend bug. **GPU kernels never compute dispersion delays themselves.** The Rust `DmPlan` computes the per-channel delay table on the host, and both the CPU and CUDA kernels consume that same table. Delay rounding is then identical by construction.

**Provenance.** Every Arrow table carries `schema_version` and provenance metadata (`software`, `git_sha`, `backend`). No astronomy project was found using Arrow for candidate exchange yet. Current practice is ad hoc: PRESTO `.singlepulse` text, heimdall `.cand` text, FETCH-style HDF5. So write export adapters to those formats rather than adopting any of them as the internal schema.

### Proposed tree

Proposed tree after the rename and restructure:

```
pulsar-transient-complete-codes/          # this repo (README.md stays)
  PULSAR-TUTORIAL/                        # submodule, unchanged
  Intro-to-Interferometry/                # submodule, unchanged
  filterbank/                             # renamed from psr-searchmode-to-filterbank
    README.md  LICENSE(MIT)  pixi.toml  Cargo.toml  justfile
    definitions/                          # THE single source of truth
      schemas/   header.yaml  dmplan.yaml  rfi_mask.yaml
                 candidate_single_pulse.yaml  candidate_periodic.yaml  fold.yaml
      codegen/   generate.py              # emits everything under generated/
      generated/ rust/  python/  cpp/include/fbdefs/*.h  arrow/*.json
    converter/                            # existing code, moved verbatim first
      pyproject.toml                      # stays independently installable
      fb_converter/ __init__.py __main__.py cli.py converter.py
                    filterbank.py psrfits.py
      scripts/fits2fil
    crates/                               # Rust CPU kernels (no Python inside)
      fb-types/        # generated structs + DM-delay math + DmPlan
      fb-io/           # fil, PSRFITS search/fold, DADA header, memmap blocks
      fb-rfi/          # IQRM, zero-DM, spectral kurtosis, normalisation
      fb-dedisp/       # brute force (reference), FDMT, subband
      fb-sp/           # prefix-sum boxcars, DBSCAN/FoF clustering
      fb-fft/          # realfft wrappers, dereddening, harmonic sum
      fb-ffa/          # FFA ported from riptide (MIT)
      fb-fold/         # batch folding, PSRFITS fold writer
      fb-py/           # the ONLY crate linking PyO3
    cuda/                                 # C++/CUDA GPU kernels
      CMakeLists.txt  include/fbcuda/*.cuh
      src/ rfi/ dedisp/ boxcar/ fdas/ harmsum/ fold/
      bindings/        # nanobind modules
      tests/           # Catch2/GoogleTest + nvbench
    python/filterbank/                    # orchestration (pure Python)
      core/            # backend registry, array contract, config
      stages/          # io, rfi, dedisp, single_pulse, periodic, ffa, fold, sift, classify, timing
      adapters/        # subprocess/optional: presto, heimdall, dmt, riptide, peasoup,
                       # pulsarx, transientx, psrchive, dspsr, pint, fetch, your
      pipelines/       # batch.py, realtime/ (ring-buffer workers), diamond.py, twoelement.py
      cli/
    realtime/          # ingest + supervision daemon (PSRDADA-compatible ring buffer)
    models/            # ONNX classifiers + training scripts
    tests/  unit/  cross_backend/  regression/  fixtures/ (git-lfs)
    bench/
    containers/        # oracle images: presto, heimdall, psrchive+dspsr, pulsarx
    docs/
```

This follows the polars and Rerun pattern: kernel crates are free of PyO3, there is one thin binding crate, and a `default-members` setting makes plain `cargo build` skip the Python bindings ([polars Cargo.toml](https://redirect.github.com/pola-rs/polars/blob/main/Cargo.toml)). Rerun also enforces dependency "bands" in CI, so a crate may depend only on crates in its own band or a lower one ([rerun ARCHITECTURE.md](https://raw.githubusercontent.com/rerun-io/rerun/main/ARCHITECTURE.md)). Copy that rule here as follows:

- `fb-types` depends on nothing.
- `fb-io` depends only on `fb-types`.
- The DSP crates depend only on `fb-types` and `fb-io`.
- `fb-py` sits on top.

Python packaging becomes three wheels:

- `filterbank`: pure Python.
- `filterbank-rs`: maturin, abi3, all operating systems.
- `filterbank-cuda-cu12` / `-cu13`: scikit-build-core, Linux only. This mirrors the RAPIDS per-CUDA-major suffix convention ([pydevtools](https://pydevtools.com/handbook/explanation/installing-cuda-python-packages/)).

Your Mac never installs anything CUDA. Check PyPI for an existing `filterbank` distribution before publishing; the import name can stay `filterbank` either way.

### Backend dispatch

Backend dispatch is a Python Protocol per stage plus a registry populated from `importlib.metadata` entry points:

```python
class Dedisperser(Protocol):
    def plan(self, hdr: FilterbankHeader, dms: DmPlan) -> Plan: ...
    def run(self, data: ArrayLike, plan: Plan, out=None) -> ArrayLike: ...

get("dedisp", backend="auto")   # cuda if the plugin is importable and a device exists, else cpu
```

The CUDA wheel registers itself only when it is installed. Each backend declares its device and dtype needs, and the orchestrator moves data only at stage boundaries. Rust kernels release the GIL so I/O can overlap with compute. CUDA entry points accept the caller's stream: nanobind accepts but does not synchronise a DLPack stream argument, and it silently copies arrays without `noconvert()` ([nanobind ndarray](https://nanobind.readthedocs.io/en/latest/ndarray.html)).

Prototype new GPU kernels as CuPy `RawKernel` strings, then promote them to compiled nanobind modules once they stabilise. **Do not route GPU work through Rust**. Doing so adds a second owner of CUDA contexts and streams. The exception is a later standalone Rust real-time binary, which can load the same `.cu`-compiled PTX through cudarc ([cudarc](https://docs.rs/crate/cudarc/latest)). CubeCL is the one "single source for CUDA and Metal" Rust option, but at 0.10/0.11 it is too young to be the primary path ([cubecl](https://docs.rs/crate/cubecl/latest)). Treat a Metal backend as optional, added kernel by kernel only where Mac speed matters.

### Real-time path

The real-time path differs from batch only in scheduling, not in the stage code. Deployed FRB systems use native ring buffers plus a control plane, not Python task schedulers. CHORD's pirate runs network ingest into a ring buffer, then RFI masking, GPU dedispersion and peak-finding, with gRPC for control ([pirate](https://kmsmith137.github.io/pirate/notes/intro.html)).

For the Diamond Array Telescope this means:

- A PSRDADA-compatible ring buffer of fixed-size gulps, overlapped by the maximum DM delay.
- Long-lived stage workers that each own a CUDA stream.
- Candidates published as Arrow RecordBatches over IPC.
- Python used only for configuration, supervision and the classifier.

Because each stage is `stage(arrays, config) -> arrays + Arrow table`, the same code runs in both batch and streaming modes. The two-element interferometer needs only an input adapter. Once its correlated or beamformed output is written as a filterbank block with a generated header, everything downstream is shared.

## Stage by stage: port the permissive leaders, call the GPL heavyweights

The table gives the verdict for each stage. "Port" means reimplement in Rust or CUDA from a permissive source or from the paper. "Wrap" means bind the library in-process (it must be permissive). "Call" means a subprocess adapter or oracle.

| Stage | Best existing implementation | Verdict | Target |
|---|---|---|---|
| I/O: fil, PSRFITS search mode | sigpyproc3 (MIT), Rust `sigproc_filterbank` (MIT/Apache) | **Port** to `fb-io` with memmap; `fitsio` for headers, direct mmap of the SUBINT rows | Rust |
| I/O: DADA, baseband (GUPPI, VDIF) | No crate exists; DSPSR is the reference | Clean-room DADA header parser; **call** DSPSR `digifil` for baseband | Rust + subprocess |
| RFI (real time) | IQRM (sigpyproc3 `iqrm_mask`), zero-DM, spectral kurtosis | **Port** all three | Rust + CUDA |
| RFI (offline benchmark) | PRESTO rfifind | **Call** as an oracle only | subprocess |
| RFI (folded data) | clfd (MIT) | **Port** the feature and Tukey logic | Rust/Python |
| Incoherent dedispersion | `dmt` FDMT and brute force (MIT) | **Wrap** first, then **port** brute force and FDMT to Rust | Python wrap, then Rust + CUDA |
| Coherent dedispersion | cdmt, DSPSR, `dmt` CohFDMT | **Call** cdmt/DSPSR; wrap `dmt` CohFDMT later | subprocess or wrap |
| Single-pulse search | TransientX design (boxcars plus DBSCAN); heimdall | **Port** from the TransientX paper (clean-room); **call** heimdall and TransientX as oracles | Rust + CUDA |
| FFT, accel and jerk search | PRESTO v6.3.3 accelsearch (CPU); PrestoZL (GPU, Fourier); peasoup (GPU, time-domain and template bank) | **Call** PRESTO and PrestoZL; **wrap or port** peasoup | subprocess; CUDA |
| FFT libraries | RustFFT/realfft (CPU), cuFFT (NVIDIA), VkFFT (Metal/portable) | Use directly | Rust / CUDA |
| FFA | riptide (MIT) | **Wrap** now, **port** the core to Rust | Python, then Rust |
| Coherent harmonic search | PRESTO `coherentsearch` (GPL), CoherentSearch.jl | **Call**; clean-room port from the paper later | subprocess |
| Folding (search candidates) | PulsarX `psrfold_fil` (GPL) | **Call** now; clean-room Rust batch folder plus CUDA later | subprocess, then Rust + CUDA |
| Folding (baseband, polarimetry) | DSPSR | **Call** | subprocess |
| Sifting | PRESTO ACCEL_sift / `coherent_sift.py` logic | Reimplement in Python over Arrow tables | Python |
| ML classification | FETCH (GPL, TF2); DRAFTS; PICS lineage | Retrain a compact PyTorch model and export ONNX; **call** FETCH as a baseline | Python + ONNX |
| Archives and timing | PSRCHIVE `pat`/`pac`, PINT | **Call** PSRCHIVE in an isolated env; **import** PINT; write your own PSRFITS fold writer | subprocess + Python |

### I/O and RFI

These are the foundation and the most bug-prone stages. No Python library covers every format: `your` reads fil and PSRFITS but can only write PSRDADA, and it is GPL ([your](https://github.com/thepetabyteproject/your)). sigpyproc3 handles fil and PSRFITS, and unpacks 1-, 2- and 4-bit data with configurable bit order, under MIT ([sigpyproc3 bits](https://sigpyproc3.readthedocs.io/en/latest/generated/sigpyproc.io.bits.html)).

In Rust, `sigproc_filterbank` is a safe, MIT/Apache parser ([docs.rs](https://docs.rs/sigproc_filterbank)). For FITS there are two choices:

- `fitsio`, which wraps cfitsio and is the most mature.
- `fitsio-pure`, which needs no C dependencies but by its own README is roughly 10x slower on large column writes ([fitsio-pure README](https://docs.rs/crate/fitsio-pure/latest/source/README.md)).

So parse headers with a crate and memory-map the fixed-width SUBINT rows directly. Sub-byte unpacking and applying DAT_SCL/DAT_OFFS/DAT_WTS are where correctness bugs cluster. Your existing `psrfits.py` and `filterbank.py` already encode this knowledge: they become the first parity test for `fb-io`, and the seed of `definitions/schemas/header.yaml`.

For RFI, use three methods:

- **IQRM.** It is non-parametric and robust to trends, and on MeerKAT and Jodrell Bank data it gave roughly 3x more genuine pulses and 3x fewer spurious candidates ([Morello RFI2022 slides](https://events.ecmwf.int/event/258/contributions/2894/attachments/1577/2816/RFI2022_Morello.pdf), [arXiv 2108.12434](https://arxiv.org/abs/2108.12434v2)). It handles narrowband interference.
- **Zero-DM filtering** ([Eatough et al. 2009](https://arxiv.org/pdf/0901.3993)) for broadband impulsive RFI. Make it switchable: it suppresses genuine low-DM signals ([arXiv 2405.20716](https://arxiv.org/pdf/2405.20716)).
- **Spectral kurtosis**, which CHIME computes on GPU for about 0.7% of X-engine compute ([arXiv 1808.10365](https://arxiv.org/pdf/1808.10365)).

Order them the way deployed systems do: channel statistics, then mask, then zero-DM, then baseline normalisation, then dedispersion. All of these are reductions over (channel, block), which map naturally to rayon and to CUDA. rfifind stays the benchmark to beat; GMRT's SPOTLIGHT reports 98% fewer false detections than rfifind ([arXiv 2608.19102](https://arxiv.org/abs/2608.19102)).

### Dedispersion and single-pulse search

`dmt` is the clear foundation:

- **Speed.** FDMT runs at **55x real time on 8 EPYC threads and 776x on one L40S**. That is about 19x (CPU) and 6x (GPU) faster than its own brute-force engine, on a 4096-channel MeerKAT setup ([pravirkr/dmt](https://github.com/pravirkr/dmt)).
- **Portability.** It has NEON SIMD, so it runs well on Apple Silicon CPUs.
- **Noise calibration.** It computes exact per-trial noise, so S/N needs no calibration run ([pravirkr/dmt](https://github.com/pravirkr/dmt)).

Wrap it in Phase 2 to get production speed immediately. Then write a Rust brute-force engine as the *reference*: it is exact and accepts an arbitrary DM grid. Barsdell et al. show the trade-offs: subband dedispersion adds about 40% extra smearing at 16 channels per subband, and about 300% at 32 ([Barsdell et al. 2012](https://ar5iv.labs.arxiv.org/html/1201.5380)). Treat brute force as truth, FDMT as the production engine, and subband (with per-range downsampling, as in TransientX) as the high-channel-count CPU fallback.

For detection, TransientX's design is proven and well documented: geometric boxcar widths, then DBSCAN with a kd-tree on the DM-time plane. It reaches real time on a single core for 4096 DMs by 4096 channels ([arXiv 2401.13834](https://arxiv.org/html/2401.13834v1)). Implement it clean-room from the paper, since the code is GPL.

The dedisp/heimdall lineage now survives only in personal forks without licence files. Heimdall++ found up to 2.66x headroom in it from engineering work alone ([arXiv 2512.00398](https://arxiv.org/abs/2512.00398)). Use heimdall as an oracle, not a foundation.

### Periodicity

The ranking of methods changed in 2026. Ransom's coherent harmonic summing recovered **76% of 860,000 injected pulsars, against 64% for riptide and 42% for accelsearch** at matched false-alarm rate ([arXiv 2610.11749](https://arxiv.org/abs/2610.11749)). PRESTO v6.2.0 ships it as `coherentsearch`, but with no acceleration dimension ([PRESTO CHANGELOG](https://github.com/scottransom/presto/blob/master/CHANGELOG.md)). PRESTO v6.3.3 also made CPU accelsearch about **8x faster** than v5.3.1 to v6.1.1 ([CHANGELOG](https://github.com/scottransom/presto/blob/master/CHANGELOG.md)). The often-quoted GPU speedups, 51 to 56x for PrestoZL and about 90x for AstroAccelerate jerk search ([PrestoZL](https://github.com/zhejianglab/PrestoZL), [arXiv 1911.01353](https://arxiv.org/pdf/1911.01353)), were measured against older PRESTO versions. Re-benchmark before investing GPU effort here.

The practical split:

- **Isolated and long-period pulsars:** call `coherentsearch`, and wrap riptide as the independent FFA baseline.
- **Binaries:** wrap peasoup on GPU. It is Apache-2.0, does time-domain resampling and template-bank search, and runs on every beam at TRAPUM ([arXiv 2212.08528](https://arxiv.org/pdf/2212.08528)). Call PrestoZL or PRESTO for Fourier-domain jerk search.

There is no Rust FFA or accelsearch anywhere, so a Rust riptide port would be new work, and the license allows it.

### Folding, classification and timing

For folding search candidates, PulsarX is about 50x faster than prepfold on large candidate sets ([arXiv 2309.02544](https://arxiv.org/pdf/2309.02544)). Call it now. Its batch, FDMT-based, regularised-fold design is also the template for a later Rust and CUDA folder.

Classifier accuracy claims need care. FETCH reports above 99.5% ([arXiv 1902.06343](https://arxiv.org/abs/1902.06343v2)), but an independent 2025 study found much lower recall on injected signals ([arXiv 2511.17797](https://arxiv.org/pdf/2511.17797)). So train your own small model on Diamond Array data, export it to ONNX so it runs on CPU, CUDA or Apple Silicon, and measure recall with your own injection harness.

PSRCHIVE's Python bindings are reported to be incompatible with NumPy 2 ([PSRCHIVE Python manual](https://psrchive.sourceforge.net/manuals/python)). Run `pat`, `pac` and `pam` as subprocesses in a separate conda environment; conda-forge now builds aarch64 ([activity](https://sourceforge.net/p/psrchive/activity)). Use PINT in-process for timing models ([arXiv 2405.01977](https://arxiv.org/html/2405.01977v2)).

## Trust comes from injection recovery against PRESTO v6.3.3 and heimdall

A new search code earns trust only by matching or beating the codes everyone already uses, on the same data. Bitwise equality is the wrong target: brute force and FDMT legitimately differ in delay quantisation. Use three test layers.

1. **Language-native unit and property tests.**
   - Rust: `cargo test` with proptest, for example "a DM-0 pulse dedisperses to itself" and "1/2/4-bit unpack round-trips".
   - CUDA: Catch2 or GoogleTest.
   - Python: pytest with hypothesis.
2. **Cross-backend equivalence.** `assert_allclose(cpu(x), cuda(x), rtol)`, parametrised over the available backends and auto-skipped without a GPU. On your Mac this runs CPU-only and stays green.
3. **Golden regression on seeded synthetic filterbanks** with injected pulses and pulsars of known DM, width, period and S/N, plus a few real files:
   - Run PRESTO (`prepsubband`, `single_pulse_search.py`, `accelsearch`, `coherentsearch`), heimdall and TransientX once, inside pinned containers.
   - Store their outputs as git-lfs fixtures.
   - Assert recall and parameter agreement within tolerances (DM bin, sample, S/N within a set percentage).

This matches what PRESTO itself now does with real Ter5 data ([PRESTO CHANGELOG](https://github.com/scottransom/presto/blob/master/CHANGELOG.md)). The Northern Cross heimdall validation sets a useful bar: its two input modes agree within 3% in TOA/DM and 10% in S/N ([arXiv 2609.20635](https://arxiv.org/abs/2609.20635)).

Benchmark with criterion (Rust), nvbench (CUDA) and pytest-benchmark (Python), with saved history. Always compare against **PRESTO v6.3.3 and `dmt`**, not against PRESTO v4.

For CI:

- Run a `macos-14` arm64 job to keep the CPU path first-class.
- Run a Linux job that compiles the CUDA code (nvcc needs no GPU) and checks that the module imports.
- Run the GPU suite nightly on a label-gated runner. GitHub-hosted GPU runners are paid, at a few cents per minute ([GitHub changelog](https://github.blog/changelog/2024-07-08-github-actions-gpu-hosted-runners-are-now-generally-available)); a self-hosted runner on your Diamond Array GPU node works too.

## Six phases from submodule rename to real-time Diamond Array triggers

### The submodule decision

Decide this first, because every later commit depends on it.

**Option A: keep `filterbank` as a submodule.** Rename the GitHub repository `katta626/psr-searchmode-to-filterbank` to `katta626/filterbank`; GitHub redirects the old URL. Then run `git mv psr-searchmode-to-filterbank filterbank` and update the URL in `.gitmodules`. All search-stack work then happens in that repository, and this repo only pins a commit.

**Option B (recommended): convert it to a regular directory.** Either use `git submodule deinit` plus `git rm`, then `git subtree add --prefix=filterbank/converter <url> <branch>`; or rewrite a clone with `git filter-repo --to-subdirectory-filter filterbank/converter` and merge it with `--allow-unrelated-histories`. Both keep the converter's history.

The case for B: shared definitions are useful only if a schema change, the code generator, the Rust crates, the CUDA headers and the Python stages can change in **one atomic commit and one CI run**. A submodule turns every such change into a two-repository dance with pinned SHAs. Your README already frames this repo as the home of your pipeline code.

The case for A: it is the right choice if other people install the converter directly from its own GitHub URL, or if you want the search stack public while the tutorials stay separate. If you choose B, archive the old repository with a pointer to the new location.

### Phases

| Phase | Scope | Concrete deliverables | Exit criterion |
|---|---|---|---|
| 0. Restructure | Submodule decision; rename to `filterbank/`; move the converter into `filterbank/converter/` unchanged; add MIT LICENSE, `pixi.toml`, `justfile`, CI skeleton (macOS arm64 + Linux) | Converter installable on its own with `pip install -e filterbank/converter`; `fits2fil` still works | `fits2fil` output on two reference PSRFITS files is byte-identical to before the move (hashes committed) |
| 1. Definitions + I/O | `definitions/schemas/*.yaml` + `codegen/generate.py`; `fb-types` (`FilterbankHeader`, `DmPlan`, delay tables); `fb-io` (fil and PSRFITS read/write, 1–32-bit, memmap blocks, DADA header); `fb-py`; `python/filterbank/core` registry | Generated Rust, Python, C++ and Arrow artefacts; CI fails on codegen drift; the converter switched to `fb-io` internally | `fb-io` matches sigpyproc3 and the old converter sample-for-sample; throughput benchmark recorded |
| 2. CPU single-pulse MVP (Mac) | `fb-rfi` (IQRM, zero-DM, SK, normalisation); `dmt` adapter plus Rust brute-force reference; `fb-sp` (boxcar + DBSCAN); Arrow candidate tables; FETCH-compatible HDF5 export; diagnostic plots; oracle containers | `filterbank search-sp file.fil` runs end to end on the Mac | Injection recovery at least matches PRESTO `single_pulse_search.py` and TransientX on the synthetic suite |
| 3. CUDA backend (Linux) | `cuda/` kernels for RFI statistics, dedispersion (wrapping `dmt` CUDA or native), boxcar; nanobind bindings; `filterbank-cuda-cu12` wheel; cross-backend tests; GPU CI | `backend="auto"` picks the GPU on the Linux host and the CPU on the Mac | CPU/GPU agree within tolerance; recall at least that of heimdall; real-time factor measured for the Diamond Array configuration |
| 4. Real-time Diamond Array | `realtime/` ingest (PSRDADA-compatible ring buffer, gulp overlap), stream workers, Arrow IPC candidate stream, ONNX classifier trained on your data, trigger and data-dump hooks, gRPC or simple control; two-element interferometer adapter | A live pipeline that keeps up with the data rate, with logged latency | Sustained above 1x real time over multi-hour runs; classifier recall measured by live injection |
| 5. Periodicity + folding | Adapters for riptide, PRESTO accelsearch/coherentsearch, PrestoZL, peasoup and PulsarX; Arrow periodic-candidate schema; Python sifting; `fb-ffa` (riptide port); `fb-fold` batch folder + PSRFITS fold writer; PSRCHIVE/PINT timing adapters | `filterbank search-periodic` produces folded, sifted, classified candidates | Recovers injected and known pulsars at least as well as the PRESTO v6.3.3 chain |
| 6. Harden + extend | cibuildwheel releases, docs, CUDA kernels for FFA and folding, optional Metal/VkFFT kernels where profiling justifies them, optional GPL add-on package if a true port is ever wanted | Tagged release | Benchmarks published against PRESTO, heimdall and `dmt` |

## Conclusion

The research shifts the question from "which giant do I port?" to "where do I draw the licence and schema boundaries?". The permissive ecosystem (`dmt`, riptide, peasoup, sigpyproc3, clfd, RustFFT, VkFFT) has matured enough to form a fast core. For example, `dmt`'s FDMT already beats brute-force GPU dedispersion by about 6x on the same card. The GPL codes stay valuable, but in a different role: oracles that prove your core is right. Two 2026 developments also warn against freezing the design too early. Coherent harmonic summing now beats both the FFA and accelsearch on sensitivity, and an 8x CPU speedup in PRESTO erased much of the old GPU advantage. A plugin registry over a stable Arrow schema is what lets you absorb the next such shift without restructuring again.

There is also a gap you are positioned to fill. No Rust dedispersion, FFA or accelsearch code exists on crates.io or GitHub. A clean, MIT-licensed Rust CPU reference, with the CPU path as the first-class backend that runs on a Mac, would be useful to the community, not only to the Diamond Array. The highest-leverage first step is still the unglamorous one: the Phase 0 move with byte-identical converter output, followed by the definitions folder that everything else will import.
