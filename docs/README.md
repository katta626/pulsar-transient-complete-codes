# filterbank docs — pulsar and transient search codes

**Start here:** [usage.md](usage.md) — build, every `fb` command, the Python API, PRESTO interoperability, CUDA.

Complete codebase for pulsar and transient searches: data conversion,
RFI mitigation, dedispersion, single-pulse and periodicity searches,
folding, candidate classification and timing.

Design: see `reference/06-our-plan.md` (the earlier `reports/unified-search-stack-plan.md` is superseded where they differ).

- **Rust** for CPU kernels, **C++/CUDA** for GPU kernels, **Python** to
  orchestrate them.
- One `definitions/` folder holds every shared type (filterbank header,
  DM plan, RFI mask, candidate schemas, fold parameters); all other
  modules depend only on it, so each stage works on its own but can
  pass data to any other.
- Every stage has a CPU backend (works on macOS) and, where useful, a
  CUDA backend (Linux with an NVIDIA GPU), chosen at run time.

## Current layout

| Path | Contents |
|---|---|
| `filterbank/converter/` | Original PSRFITS → filterbank converter (`fits2fil`), see `converter.md` |
| `filterbank/crates/` | Rust: `fb-defs` (shared definitions), `fb-io` (PSRFITS/SIGPROC), `fb-formats` (PRESTO, heimdall, candidate files) |
| `docs/reference/` | Reference set 00–06: overview, formats, PRESTO, DSPSR/PSRCHIVE/timing, other codes, math, our plan |
| `docs/reports/`, `docs/research_notes/` | Earlier report and the source notes behind the reference set |
| [`docs/CREDITS.md`](CREDITS.md), [`docs/math-maps/`](math-maps/README.md) | What our code took from where; one math map per source library |

## Roadmap

0. Restructure the repo, keep converter output byte-identical — done.
1. Shared definitions and Rust I/O (PSRFITS, SIGPROC, DADA, PRESTO formats) — done.
2. CPU single-pulse pipeline (RFI, brute/subband/FDMT dedispersion, boxcar, clustering) — done.
3. CUDA backend — kernels for dedispersion, FDMT, zero-DM, boxcar and FFT periodicity search. Every kernel matches the CPU code under the CPU emulator (`tests/test_cuda_emu.py`); not yet compiled with `nvcc` or run on a GPU.
4. Real-time Diamond Array pipeline — growing-file follower, PSRDADA shared-memory ring reader/writer and classifier done.
5. Periodicity search and folding (FFT, time-domain and Fourier-domain acceleration, jerk search with candidate optimisation, FFA, sift, fold, TOAs) — done.
6. Vendoring of upstream engines (dmt, PRESTO accelsearch/coherentsearch, riptide core), benchmarks against the originals.
