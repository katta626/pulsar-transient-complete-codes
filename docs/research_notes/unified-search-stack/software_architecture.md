# Software Architecture for a Polyglot (Rust / C++-CUDA / Python) Pulsar and FRB Search Stack

Scope note: research done 2026-11-24. Items under "Cited Findings" are backed by a fetched or searched source. Items under "Inferences" are the researcher's engineering judgement, partly from background knowledge that was not re-checked in this session. Treat them as recommendations, not verified facts.

## Q1. Python <-> Rust bindings (PyO3, maturin, rust-numpy, zero-copy)

### Takeaway
PyO3 (0.29.x line as of late 2026) plus maturin is the standard path. Free-threaded CPython is supported, but abi3 wheels cannot load on free-threaded builds, so plan per-version wheels for 3.14t. For arrays, use rust-numpy for NumPy and pyo3-arrow (Arrow PyCapsule interface) for tables such as candidate lists. Both give zero-copy at the buffer level.

### Cited Findings
- PyO3's current docs are on the 0.29.x line (0.29.0 to 0.29.3 on docs.rs). — [docs.rs pyo3 0.29.3](https://docs.rs/crate/pyo3/0.29.3)
- Free-threaded builds cannot load `abi3` wheels. Both GIL and free-threaded builds can load `abi3t` wheels, but `abi3t` targets only Python 3.15 and newer. On 3.14 an `abi3t` request quietly produces a version-specific extension. — [PyO3 guide: building and distribution, 0.29.1](https://docs.rs/crate/pyo3/0.29.1/source/guide/src/building-and-distribution.md); [PyO3 features 0.29.0](https://docs.rs/crate/pyo3/0.29.0/source/guide/src/features.md)
- A downstream maintainer confirms the practical effect: 3.14 free-threaded support needs version-specific wheels. — [openai/harmony issue 87](https://github.com/openai/harmony/issues/87)
- Free-threaded support first landed around PyO3 0.23. Since 0.28, PyO3 assumes modules are thread-safe by default, and a module can opt out with `gil_used = true`. — [PyO3 free-threading guide 0.23](https://docs.rs/crate/pyo3/0.23.0/source/guide/src/free-threading.md); [0.28](https://docs.rs/crate/pyo3/0.28.0/source/guide/src/free-threading.md)
- pyo3-arrow (0.19.0 latest) does zero-copy FFI between Python Arrow objects and Rust `arrow` objects through the Arrow PyCapsule Interface. Its structs used as function arguments auto-convert input. Buffer-protocol (numpy) input is also zero-copy. — [docs.rs pyo3-arrow](https://docs.rs/crate/pyo3-arrow/latest)
- The Arrow PyCapsule Interface standardizes ArrowSchema, ArrowArray, and ArrowArrayStream capsules. It is marked experimental, and a capsule can only be consumed once. — [Apache Arrow PyCapsule Interface](https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html)
- Caveat: some crates avoid pyo3-arrow because of compile or version friction (it has had many versions, 0.6 to 0.19) and implement the PyCapsule protocol directly on arrow-array plus pyo3. Utf8View import is a copy, not zero-copy. — [minarrow-pyo3 README](https://docs.rs/crate/minarrow-pyo3/latest/source/README.md)
- Polars is the reference layout. The Python bindings are a separate crate (`py-polars`, now `polars-python`, PyO3 with abi3) that depends on the core `polars-*` crates through workspace references. It is built with `cd py-polars && maturin develop --release`. — [polars-python Cargo.toml](https://docs.rs/crate/polars-python/0.54.4/source/Cargo.toml.orig); [polars Cargo.toml](https://docs.rs/crate/polars/latest/source/Cargo.toml.orig)
- DLPack in Rust: the `dlpark` crate (0.8, 0.9 alpha) provides DLPack import and export for PyO3, which can carry CPU or CUDA tensors. — [dlpark](https://docs.rs/crate/dlpark/0.8.0)

### Inferences
- Pattern: keep pure-Rust kernel crates (`psr-dedisp`, `psr-fft`, etc.) free of PyO3, and add one thin `*-py` binding crate per Python package. This is the polars pattern, and it lets Rust users and tests avoid linking libpython.
- Use rust-numpy `PyReadonlyArray2<f32>` / `PyReadwriteArray` for borrowed views, and release the GIL (`py.allow_threads` / `Python::detach` in newer PyO3) around long kernels so the Python orchestrator can overlap I/O.
- Ship abi3 wheels for 3.10 to 3.14 GIL builds (one wheel per platform), plus version-specific `cp314t` wheels if free-threading matters. Revisit `abi3t` once 3.15 is the floor.
- The rust-numpy version must match PyO3 (rust-numpy tracks PyO3 minor versions). Pin both together in the workspace `[workspace.dependencies]`.

### Gaps
- rust-numpy's current version and its behaviour under free-threaded CPython (borrow-checking of arrays across threads) were not verified in this session.

## Q2. Python <-> C++/CUDA bindings (nanobind, pybind11, scikit-build-core, CuPy, DLPack, CAI)

### Takeaway
For new C++/CUDA code, use nanobind built by scikit-build-core. `nb::ndarray` gives typed, device-constrained zero-copy views over NumPy, CuPy, PyTorch, JAX, and MLX arrays through DLPack. CuPy RawKernel/RawModule is a lighter path for prototyping kernels with no compiled extension at all.

### Cited Findings
- `nb::ndarray` exchanges arrays zero-copy with NumPy, PyTorch, TensorFlow, JAX, CuPy, and MLX through the buffer protocol or DLPack. Template parameters constrain dtype, shape, device (`nb::device::cpu`, `nb::device::cuda`), and order (`nb::c_contig`). Return-type tags include `nb::numpy`, `nb::cupy`, `nb::pytorch`, `nb::jax`, `nb::mlx`, and `nb::array_api` (which exposes `__dlpack__` and `__dlpack_device__`). — [nanobind ndarray docs](https://nanobind.readthedocs.io/en/latest/ndarray.html)
- Pitfalls in the nanobind docs:
  - Returning an array without an `owner` or capsule makes nanobind assume the memory lives forever, which is a common use-after-free.
  - nanobind accepts the `stream` argument of `__dlpack__` but does not synchronize it. Asynchronous device output needs explicit stream sync.
  - Without `noconvert()`, nanobind may silently copy to fix dtype or layout.
  — [nanobind ndarray docs](https://nanobind.readthedocs.io/en/latest/ndarray.html)
- nanobind's packaging guide uses scikit-build-core (>=0.4.3) plus nanobind in `pyproject.toml` and recommends stable-ABI wheels for CPython 3.12+ via `wheel.py-api = "cp312"`. — [nanobind packaging](https://nanobind.readthedocs.io/en/latest/packaging.html)
- Recent nanobind release line is 2.9.x (Debian ships 2.9.2). — [Debian nanobind source](https://sources.debian.org/src/nanobind/2.9.2-2/docs/ndarray.rst)
- NVIDIA's own `cuda.core` (cuda-python) has a DLPack and memory-view layer (`_dlpack.pyx`, `_memoryview.pyx`) and is installed as `pip install cuda-core[cu12]`. It is now past 1.0 (docs for 1.0.0 and 1.2.1 exist). — [cuda.core 0.7.0 notes, 1.2.1 docs tree](https://nvidia.github.io/cuda-python/cuda-core/1.2.1/release/0.7.0-notes.html); [cuda-core on PyPI](https://pypi.org/project/cuda-core/)
- Domain precedent: the CHORD real-time FRB search "pirate" is C++/CUDA with pybind11 bindings (`pirate_frb`). Complex GPU kernels are generated by a Python code generator, and the Python surface is limited to high-level configuration. It builds only on Ubuntu with a physical GPU. — [pirate intro notes](https://kmsmith137.github.io/pirate/notes/intro.html)

### Inferences
- nanobind vs pybind11: nanobind has smaller binaries, faster compiles, lower call overhead, and native DLPack ndarray. pybind11 is still fine when wrapping existing pybind11 code (pirate, many astro codes). For new wrappers of heimdall or dedisp-style kernels, prefer nanobind.
- CuPy `RawKernel` / `RawModule` (NVRTC at runtime) is ideal for prototype GPU kernels: no build step, works from a `.cu` string, and interoperates through `__cuda_array_interface__` and DLPack. Promote a kernel to a compiled nanobind module once it stabilizes, or when you need multi-kernel C++ host logic, cuFFT plans, or CUDA graphs.
- Interop rule: every kernel entry point should accept any DLPack or CAI producer and return a DLPack-capable object. Then NumPy, CuPy, PyTorch, and JAX users all work without copies. Always honour or perform stream synchronization at module boundaries.
- In scikit-build-core, use `CMAKE_CUDA_ARCHITECTURES` with a fat binary (for example `75;80;86;89;90;100;120` plus PTX for the newest arch) and `wheel.py-api = "cp312"` for stable-ABI wheels.

### Gaps
- The exact DLPack version (1.0 / 1.1, including the versioned `__dlpack__(max_version=...)`) targeted by current nanobind was not stated on the fetched page.
- The current CuPy version and its RawKernel API were not re-verified this session.

## Q3. Rust <-> CUDA (cudarc, cust, rust-cuda), and routing GPU through Rust or directly from Python

### Takeaway
cudarc is the maintained, pragmatic Rust CUDA crate (0.19.x): it covers the driver, NVRTC, cuBLAS, cuRAND, and more, with sys / result / safe layers. Portable-kernel Rust (CubeCL) is promising but young. Recommendation: write CUDA kernels in .cu (C++), expose them to Python with nanobind, and expose them to Rust only when a Rust host pipeline really needs them (through cudarc loading PTX or cubin, or a C ABI). Do not force all GPU work through Rust.

### Cited Findings
- cudarc 0.19.9 is the latest. It provides safe wrappers for the CUDA driver, NVRTC, cuRAND, cuBLAS, cuBLASLt, and others. `nvrtc::compile_ptx()` produces PTX that is loaded with `CudaContext::load_module()`. Each wrapper has three layers: sys (bindgen FFI), result, and safe. — [cudarc docs.rs](https://docs.rs/crate/cudarc/latest); [README](https://docs.rs/crate/cudarc/0.19.6/source/README.md)
- The CUDA version is chosen with Cargo features (for example `cuda-12010`), or detected with `cuda-version-from-build-system`. — [cudarc Cargo features](https://docs.rs/crate/cudarc/0.19.1/source/Cargo.toml)
- The repo moved to `github.com/chelsea0x3b/cudarc`. — [GitHub](https://github.com/chelsea0x3b/cudarc)
- CubeCL (0.10.0 stable, 0.11 pre-releases) writes kernels as Rust functions marked `#[cube]` and JIT-compiles them to backends:
  - WebGPU/wgpu (WGSL)
  - CUDA
  - ROCm/HIP
  - Metal (C++ Metal compiler, on Apple GPUs)
  - Vulkan
  - CPU

  It offers autotuning and comptime specialization. Feature support varies by backend, for example Tensor Cores are not available on WebGPU. — [cubecl docs.rs](https://docs.rs/crate/cubecl/latest); [tracel-ai/cubecl](https://github.com/tracel-ai/cubecl)

### Inferences
- `cust` and Rust-CUDA (`rustc_codegen_nvvm`): Rust-CUDA was revived around 2025, but it still needs a pinned nightly toolchain and has a small user base. Not recommended as the primary path for a science codebase. (Not re-verified this session.)
- Routing decision:
  - Python is the orchestrator and owns the device buffers (CuPy or PyTorch arrays).
  - GPU kernels are C++/CUDA modules reached directly from Python through nanobind and DLPack.
  - Rust is for CPU kernels and I/O (filterbank and PSRFITS readers, CPU dedispersion, FFA, sifting).
  - Routing GPU through Rust adds a second FFI hop and a second owner of CUDA contexts and streams, which risks context or stream mismatch with CuPy or PyTorch. It is only worth it if a standalone Rust binary (no Python) must run on GPU, for example a real-time daemon.
- If a Rust-native GPU path is wanted later, the cleanest option is for the same `.cu` sources to emit PTX or fatbin at build time, with Rust loading them through cudarc. That keeps one kernel source.

### Gaps
- No benchmark comparing CubeCL with hand-written CUDA for dedispersion or FFT-heavy workloads was found.

## Q4. Shared "definitions" layer (single source of truth for types)

### Takeaway
Use a two-tier approach:
- Bulk data (dynamic spectra, time series): plain N-d arrays exchanged through DLPack or the buffer protocol.
- Records (filterbank header, DM plan, candidate table): an Arrow schema, so a candidate list is an Arrow RecordBatch readable zero-copy in Rust (`arrow`), Python (pyarrow, polars), and C++ (Arrow C Data Interface).

Generate small struct types (header, plan) from one IDL. Rerun's model is the best precedent: one definition crate generates Rust, Python, and C++.

### Cited Findings
- Rerun defines its types in a Rust-subset IDL (`crates/build/re_type_definitions`, `#[rerun_type]` macro). A code generator (`crates/build/re_types_builder`) emits the Rust, Python, and C++ SDKs. A separate `re_protos_builder` generates Protobuf and gRPC code. All log data, network transport, `.rrd` files, and the in-memory store use Apache Arrow (`re_chunk`). — [rerun ARCHITECTURE.md](https://raw.githubusercontent.com/rerun-io/rerun/main/ARCHITECTURE.md)
- Pirate defines its upstream metadata in YAML (`configs/xengine_metadata.yml`), uses gRPC for control and monitoring, and writes data with asdf-cxx. — [pirate intro](https://kmsmith137.github.io/pirate/notes/intro.html)
- The Arrow PyCapsule Interface lets any Arrow-speaking library (pyarrow, polars, arro3, Rust arrow through pyo3-arrow) exchange tables without copies. — [Arrow PyCapsule Interface](https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html); [pyo3-arrow](https://docs.rs/crate/pyo3-arrow/latest)
- An existing pure-Rust SIGPROC filterbank crate exists: `sigproc_filterbank` (kiranshila), no unsafe, read/write only. — [sigproc_filterbank](https://docs.rs/crate/sigproc_filterbank/0.3.0)

### Inferences
Options compared:

| Option | Strengths | Weaknesses | Verdict |
|---|---|---|---|
| Arrow schema (records and tables) | Zero-copy everywhere, columnar, Parquet and IPC on disk, native polars and pandas | Not great for small nested config | Use for candidates, DM-trial metadata, sifted lists |
| Rust core crate + PyO3 classes + cbindgen C header | One source in Rust, real methods (validation, DM-delay math) | C++ consumes a C ABI, not idiomatic C++; Python types are opaque unless also exposed as dataclasses | Good for the header and plan structs plus logic |
| FlatBuffers / Cap'n Proto | Zero-copy wire format, codegen for all 3 languages | Adds IDL toolchain; awkward Python ergonomics; overkill unless streaming over the network | Only for a real-time network protocol |
| Protobuf + gRPC | Mature, good for control-plane RPC | Copies on parse; not for bulk data | Use for real-time control and monitoring (as pirate does) |
| JSON Schema / YAML + codegen (pydantic + serde + C++ struct) | Human-editable configs | Codegen glue to maintain | Good for pipeline config files |

- Recommended concrete design:
  - A `psr-types` Rust crate holds `FilterbankHeader`, `DmPlan`, `Candidate`, and enums for units and dtypes, with `serde`, plus an `arrow_schema()` function for each tabular type.
  - A small codegen script emits `candidate.arrow.json` (the schema), Python `msgspec` or `pydantic` dataclasses, and a C++ header (cbindgen for `#[repr(C)]` POD structs).
  - CI asserts that the generated files are up to date (a "codegen diff must be empty" check, as Rerun's CI enforces layering).
- Keep the candidate schema columnar and versioned. Include a `schema_version` key in the Arrow metadata, plus provenance metadata (`software`, `git_sha`, `backend=cpu|cuda`).

### Gaps
- No astronomy project was found that uses Arrow for single-pulse candidate exchange. Practice is still ad hoc: PRESTO `.singlepulse` text, heimdall `.cand` text, and HDF5 or FETCH-style files.

## Q5. Monorepo layout patterns (polars, rerun, pirate, RAPIDS)

### Takeaway
Use a Cargo workspace under `crates/` with layered bands, separate binding crates, a CMake/scikit-build-core subtree for CUDA, and Python packages under `python/`. A single task runner (pixi or uv plus just) drives all three toolchains.

### Cited Findings
- Polars: the workspace members are `crates/*`, `docs/src/rust`, and `py-polars`. `default-members` is only `crates/*`, so `cargo build` skips the Python bindings. Resolver is v2. — [polars Cargo.toml](https://redirect.github.com/pola-rs/polars/blob/main/Cargo.toml)
- Rerun:
  - Crates are grouped in bands (folders) under `crates/`. A crate may depend only on crates in its own band or a lower one. `scripts/check_crate_dependencies.py` enforces this in CI.
  - `rerun_py` and `rerun_c` are in the top band.
  - — [rerun ARCHITECTURE.md](https://raw.githubusercontent.com/rerun-io/rerun/main/ARCHITECTURE.md)
- Rerun uses pixi tasks (`py-build`, `py-test`, `py-lint`) as the single dev entry point. pixi generates `rerun_py/pyo3-build.cfg` so Cargo caching stays consistent across build methods. — [rerun BUILD.md](https://raw.githubusercontent.com/rerun-io/rerun/main/BUILD.md)
- RAPIDS publishes per-CUDA-major packages (`cudf-cu12`, `cudf-cu13`). cuML wheels have been on PyPI since 25.10. — [pydevtools CUDA packages](https://pydevtools.com/handbook/explanation/installing-cuda-python-packages/); [NVIDIA cuML on PyPI blog](https://developer.nvidia.com/blog/reducing-cuda-binary-size-to-distribute-cuml-on-pypi)

### Inferences
Recommended 2026 layout:
```
psrstack/
  Cargo.toml                 # [workspace] members=crates/*/*, default-members exclude *-py
  pixi.toml                  # conda-forge env: rust, cmake, ninja, cuda-toolkit (linux), python; tasks
  crates/
    core/   psr-types, psr-units, psr-plan            # band 0: no deps on others
    io/     psr-filterbank, psr-psrfits               # band 1
    dsp/    psr-dedisp-cpu, psr-ffa, psr-fft, psr-rfi # band 2 (rayon, rustfft, ndarray)
    py/     psr-io-py, psr-dsp-py                     # band 3: PyO3 bindings only
  cuda/
    CMakeLists.txt           # project(LANGUAGES CXX CUDA); one target per kernel family
    include/psr_cuda/*.cuh   # C++ API, no Python
    src/dedisp/*.cu  src/boxcar/*.cu  src/fdas/*.cu
    bindings/*.cpp           # nanobind modules
    tests/ (Catch2 or GoogleTest)
  python/
    psrstack/                # pure-python orchestrator, backend registry, CLI
      pyproject.toml (hatchling)
    psrstack-rs/             # maturin; manifest-path = ../../crates/py/psr-dsp-py
      pyproject.toml
    psrstack-cuda/           # scikit-build-core; cmake.source-dir = ../../cuda
      pyproject.toml
  schemas/                   # generated arrow schemas + codegen script
  third_party/               # git submodules: presto, heimdall (reference only, GPL-isolated)
  tests/regression/          # golden outputs vs PRESTO/heimdall
  bench/
```
- Three wheels: `psrstack` (pure Python, depends on `psrstack-rs`), `psrstack-rs` (abi3, all OSes), and `psrstack-cuda-cu12` / `-cu13` (Linux x86_64 and aarch64 only, optional extra `psrstack[cuda12]`). This mirrors the RAPIDS suffix convention and `cuda-core[cu12]` extras. A Mac user never installs anything CUDA.
- Use a uv workspace or pixi for dev. pixi is preferable here because conda-forge provides `cuda-toolkit`, `cuda-nvcc`, and compilers reproducibly on Linux, and plain clang on macOS.

### Gaps
- The CuPy and JAX repo layouts were not fetched in this session. They use setup and Bazel-style builds that are less relevant.

## Q6. GPU portability with macOS Apple Silicon development

### Takeaway
Develop and unit-test on Mac against the CPU backend, which is the reference implementation. Run CUDA on Linux servers and CI. Abstract backends behind a Python-level registry with a stable array contract (DLPack). A Metal or wgpu backend is optional and possible via CubeCL or MLX, but not worth maintaining for every kernel.

### Cited Findings
- CubeCL targets CUDA, HIP, Metal, Vulkan, WebGPU, and CPU from one Rust kernel source. This is the strongest "one source, Mac and NVIDIA" option in Rust, but it is young (0.10 / 0.11-pre). — [cubecl](https://docs.rs/crate/cubecl/latest)
- nanobind can return `nb::mlx` arrays, so a C++ module can hand data to Apple MLX on Mac. — [nanobind ndarray](https://nanobind.readthedocs.io/en/latest/ndarray.html)
- uv's `--torch-backend=auto` detects visible drivers and falls back to CPU. OS markers cannot detect GPUs, so Linux CPU hosts receive CUDA wheels under a naive mapping. — [pydevtools](https://pydevtools.com/handbook/explanation/installing-cuda-python-packages/)

### Inferences
- Options ranked for this project:
  1. CUDA, plus Rust/rayon CPU fallback. This is the recommended default.
  2. CuPy with a NumPy fallback via the array-API namespace (`array_api_compat.array_namespace(x)`) for glue code.
  3. CubeCL for new portable kernels where Metal speed matters.
  4. HIP via hipify, only if AMD clusters appear.
  5. SYCL and Kokkos are heavy and C++-centric. Avoid them unless joining an HPC centre standard.
  6. Raw Metal or wgpu. Avoid, since it is per-kernel duplication.
- Backend-dispatch interface (Python):
  ```python
  class Dedisperser(Protocol):
      def plan(self, hdr: FilterbankHeader, dms: DmPlan) -> Plan: ...
      def run(self, data: ArrayLike, plan: Plan, out=None) -> ArrayLike: ...
  BACKENDS = {"cpu": RustDedisperser, "cuda": CudaDedisperser}  # entry points
  def get(kind, backend="auto"): ...  # auto = cuda if importable & device present else cpu
  ```
  - Register backends through `importlib.metadata` entry points (`psrstack.backends`), so `psrstack-cuda` plugs in only when installed.
  - Every backend takes and returns DLPack-capable arrays. The backend reports its `device` and `dtype` requirements, and the orchestrator moves data only at stage boundaries.
- Keep numerics identical across backends: same DM-delay rounding, same float32 accumulation order where feasible. Define a tolerance per stage in the regression tests.

### Gaps
- No measured Metal-vs-CUDA dedispersion performance for CubeCL or MLX was found.

## Q7. Pipeline orchestration for streaming and real-time

### Takeaway
Real-time FRB systems such as pirate and heimdall-based ones use native ring buffers plus a control plane (gRPC), not Python task schedulers. Batch archival searches can use a simple DAG, or Dask or Ray for fan-out over files and DM ranges.

### Cited Findings
- Pirate:
  - Network ingest (TCP from the X-engine) feeds a ring buffer, then RFI masking, GPU dedispersion, and peak-finding normalized to SNR.
  - Disk writes go in two stages: local SSD, then NFS.
  - gRPC services handle control, monitoring, and triggered dumps, with a Python client.
  - — [pirate intro](https://kmsmith137.github.io/pirate/notes/intro.html)
- Heimdall++ reports up to 2.66x single-file and 2.05x multi-file speedups over Heimdall by improving GPU utilization and pipeline parallelism, with results consistent with the original. — [Heimdall++ arXiv 2512.00398](https://arxiv.org/html/2512.00398v1)
- ASTROFLOW is a 2025 real-time end-to-end single-pulse pipeline (paper exists, architecture not fetched). — [arXiv 2511.02328](https://arxiv.org/html/2511.02328v1)

### Inferences
- Real-time path:
  - Ingest goes to a lock-free SPSC ring buffer of fixed-size gulps (shared memory, or PSRDADA for compatibility with existing backends).
  - Stages are long-lived workers that each own a CUDA stream. Gulps overlap by the maximum DM delay.
  - The candidate stream is Arrow RecordBatches over IPC.
  - Python is used only for config and supervision.
- Batch path: a plain Python DAG of pure functions (file to header to plan to dedisp to search to sift to fold) is enough. Wrap it with Prefect or Dask only for cluster fan-out and retries. Ray suits heterogeneous CPU/GPU actor pools.
- Design every stage as `stage(input_arrays, config) -> output_arrays + Arrow table`. Then the same stage runs in batch or inside a streaming worker loop.

### Gaps
- Fresh 2026 comparisons of Prefect, Dask, and Ray for HPC/Slurm astronomy pipelines were not researched.

## Q8. Testing and benchmarking

### Takeaway
Use a three-layer test pyramid: language-native unit tests, cross-backend equivalence tests, and golden regression against PRESTO or heimdall outputs on fixed synthetic and real datasets. GPU tests run on a label-gated GPU runner, everything else on hosted CPU runners.

### Cited Findings
- GitHub-hosted GPU runners (T4) have been GA since 2024 for Linux and Windows. They are paid, at a few cents per minute, with no free tier. — [GitHub changelog](https://github.blog/changelog/2024-07-08-github-actions-gpu-hosted-runners-are-now-generally-available); [scientific-python blog](https://blog.scientific-python.org/scikit-learn/gpu-ci/)
- scikit-learn triggers GPU CI with a "CUDA CI" PR label. RAPIDS uses NVIDIA-provided self-hosted runners in a `linux-nvidia-gpu` group with label gating. — [scientific-python blog](https://blog.scientific-python.org/scikit-learn/gpu-ci/); [NVIDIA GitHub Actions docs](https://docs.nvidia.com/datascience/deployment/latest/developer/ci/github-actions/)
- The common pattern is to build CUDA wheels with cibuildwheel in manylinux images that include the CUDA toolkit, import-check on CPU runners, and run kernels only on GPU machines. — [cibuildwheel](https://github.com/pypa/cibuildwheel); [manylinux-cuda](https://github.com/ameli/manylinux-cuda)
- Alternatives include Depot GPU runners and Hugging Face Jobs as GitHub Actions executors. — [Depot](https://depot.dev/blog/github-actions-gpus-now-available); [HF blog](https://huggingface.co/blog/github-ci-hf-jobs)

### Inferences
- Rust: `cargo test` plus `proptest` for invariants (dedispersion of a zero-DM pulse is the identity), and `criterion` for microbenchmarks. C++: Catch2 or GoogleTest plus nvbench or Google Benchmark. Python: pytest, pytest-benchmark, and `hypothesis`.
- Golden regression:
  - Generate deterministic synthetic filterbanks with injected pulses and pulsars of known DM, width, and SNR, for example with a seeded injector.
  - Run PRESTO (`prepsubband`, `single_pulse_search.py`, `accelsearch`) and heimdall once in a container, and commit their outputs or hashes as fixtures. Use git-lfs or a data bucket.
  - Assert recall and matching within tolerances (DM bin, time sample, SNR within x%), not bitwise equality.
- Cross-backend test: `assert_allclose(cpu(x), cuda(x), rtol)` parametrized over the backends available, auto-skipped when CUDA is absent. On Mac this silently runs CPU only.
- Store benchmark history with pytest-benchmark `--benchmark-autosave`, or Bencher or codspeed for regression alerts.

### Gaps
- No public benchmark suite for pulsar search codes exists to reuse. Synthetic data generators (simulateSearch, injectpulsar-style tools) were not surveyed here.

## Q9. Licensing implications of porting GPL / AFL code

### Takeaway
PRESTO is GPL-2.0. A line-by-line port or translation into Rust is very likely a derivative work and must stay GPL-compatible. DSPSR and PSRCHIVE are listed under the Academic Free License, and heimdall's licence was not verified. Either license the whole stack GPL, or do clean-room reimplementations from papers with PRESTO used only as an external test oracle.

### Cited Findings
- PRESTO's licence is GPL-2.0 (COPYING file). — [scottransom/presto](https://github.com/scottransom/presto); [ecosyste.ms](https://awesome.ecosyste.ms/projects/github.com%2Fscottransom%2Fpresto)
- PSRCHIVE and DSPSR are listed on SourceForge as Academic Free License (AFL). — [psrchive SourceForge](https://sourceforge.net/projects/psrchive/); [dspsr SourceForge](https://sourceforge.net/projects/dspsr/)
- The FSF GPL FAQ has an entry titled "What does the GPL say about translating some code to a different programming language?". A secondhand 2018 summary says it treats translation as requiring the program's licence permission, that is, a translation is a modified version. The FAQ page could not be fetched (HTTP 429), so the exact wording was not confirmed. — [GNU GPL FAQ](https://www.gnu.org/licenses/gpl-faq.html); [Slashdot summary](https://news.slashdot.org/story/18/11/09/2311224/the-free-software-foundation-releases-new-comments-about-licenses)
- Wrapper projects around heimdall carry their own licences (for example hdpipe is MIT), which covers only the wrapper. — [science.ecosyste.ms entry](https://science.ecosyste.ms/projects/51633)

### Inferences (not legal advice)
- If you port PRESTO C to Rust, the result is GPLv2 (PRESTO lacks an "or later" clause, so check the headers). Linking it into a PyO3 module makes that wheel GPL. Keep GPL-derived code in separately distributed crates or wheels (`psrstack-presto-port`) so the MIT/Apache core stays permissive, or license everything GPL-2.0-or-later-compatible from the start.
- Using PRESTO or heimdall as subprocess test oracles, or reading their file formats, does not make your code derivative.
- Clean-room reimplementation from published algorithms (Ransom 2002 accel search, Barsdell 2012 dedispersion, FFA papers) avoids GPL inheritance. Document the provenance of each module.
- AFL-3.0 is OSI-approved but is generally regarded as GPL-incompatible. Do not mix AFL-derived (DSPSR/PSRCHIVE) code and GPL-derived (PRESTO) code in the same binary without checking.
- Apache-2.0 is incompatible with GPLv2-only, while MIT is compatible. If the core is permissive and might ever link GPLv2 parts, prefer MIT or BSD.

### Gaps
- Heimdall's licence was not found. Check `LICENSE` in the heimdall-astro SourceForge or GitHub repo.
- The exact FSF FAQ wording and whether PRESTO files say "version 2 or later" were not verified.

## Cross-cutting key questions

### Takeaway
Independent-yet-interoperable modules come from four things:
1. Each module is a function over arrays plus Arrow tables, with no hidden global state.
2. There is one shared types crate or schema.
3. Every language boundary is crossed through zero-copy standard protocols (DLPack/CAI for arrays, Arrow C Data/PyCapsule for tables).
4. Backends are optional plugins discovered at runtime.

The biggest pitfalls are CUDA wheel distribution, ABI (abi3 vs free-threaded), and GPU-less CI.

### Cited Findings
- PEP 817 (wheel-variant model) is a draft and PEP 825 (variant format) is provisionally accepted. Installers do not yet auto-select GPU variants by default; Astral has an experimental uv build. Today's practice is suffixed names (`-cu12`), extras (`jax[cuda13]`, `cuda-core[cu12]`), or separate indexes (PyTorch). — [pydevtools](https://pydevtools.com/handbook/explanation/installing-cuda-python-packages/); [PEP 817](https://peps.python.org/pep-0817/); [PEP 825](https://peps.python.org/pep-0825/); [NVIDIA wheel variants blog](https://developer.nvidia.com/blog/streamline-cuda-accelerated-python-install-and-packaging-workflows-with-wheel-variants); [Astral wheel variants](https://astral.sh/blog/wheel-variants)
- CUDA wheel pitfalls:
  - Wheels are huge, with a 1 GB hard cap and PyPI project limits.
  - Wheels lack symlinks, so `nvidia-*` runtime wheels need special handling.
  - pip cannot express a dependency on the CUDA driver or detect GPUs.
  - Build per-CUDA-major wheels to use minor-version compatibility.
  - PTX JIT is not covered by minor-version compatibility, so PTX from a newer toolkit can fail on older drivers; nvJitLink mitigates this.
  - Shipping CUDA from conda is an established alternative.
  - — [pypackaging-native: GPUs](https://pypackaging-native.github.io/key-issues/gpus)

### Inferences
- Pitfall checklist:
  - **Duplicate CUDA runtimes or contexts.** Link the CUDA runtime statically (`-cudart static`), or depend on `nvidia-cuda-runtime-cu12` wheels. Share the primary context and use the caller's stream: accept a stream argument and default to CuPy or PyTorch's current stream.
  - **Fat binaries.** Compile cubins for the supported compute capabilities plus PTX for the newest one only, to keep wheels small. Consider loading cuFFT from `nvidia-cufft-cu12` instead of bundling it.
  - **manylinux.** Use manylinux_2_28 images with CUDA preinstalled, run `auditwheel repair --exclude libcuda.so.1` (the driver must never be bundled), and on macOS ship no CUDA wheel.
  - **Rust and C++ in one process.** Avoid two copies of the same C library (for example FFTW or cfitsio) statically linked into different extension modules. Prefer one I/O implementation, in Rust.
  - **abi3.** nanobind uses `cp312` stable ABI and PyO3 uses `abi3-py310`. Free-threaded 3.14t needs separate wheels for both.
  - **CI without GPUs.** On hosted runners, build CUDA code (nvcc compiles without a GPU), import the module, and skip device tests via `pytest.mark.cuda` and a `has_cuda()` fixture. On label-gated self-hosted or paid GPU runners, run the full suite nightly and on release.
  - **Mac development.** `pixi run test` must pass with zero CUDA. CI includes a `macos-14` arm64 job to keep the CPU path first-class.
- Build tooling for 2026:
  - Cargo workspace and maturin (>=1.x) for Rust wheels.
  - CMake (>=3.28) with scikit-build-core and nanobind for CUDA wheels.
  - pixi (conda-forge) or uv for environments and tasks.
  - cibuildwheel for release builds.
  - ruff and mypy/pyright for Python, clippy and rustfmt for Rust, clang-format for C++.
  - pre-commit to enforce generated-schema freshness.

### Gaps
- No existing open-source astronomy project combining Rust, C++/CUDA, and Python in one monorepo was found. The closest precedents are rerun (Rust, C++, Python, non-astronomy), polars (Rust and Python), and pirate (C++/CUDA and Python, FRB). Rust pulsar crates found are small I/O utilities (`sigproc_filterbank`, `psrutils`), with no Rust dedispersion or search crates. — [sigproc_filterbank](https://docs.rs/crate/sigproc_filterbank/0.3.0); [psrutils](https://docs.rs/crate/psrutils/0.2.2/source/Cargo.toml.orig)
