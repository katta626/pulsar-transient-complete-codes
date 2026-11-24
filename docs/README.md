# filterbank docs — pulsar and transient search codes

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

## Roadmap

0. Restructure the repo, keep converter output byte-identical (done).
1. Shared definitions and Rust I/O.
2. CPU single-pulse pipeline.
3. CUDA backend.
4. Real-time Diamond Array pipeline.
5. Periodicity search and folding.
6. Release hardening.
