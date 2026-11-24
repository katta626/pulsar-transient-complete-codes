# Pulsar/Transient Data I/O Formats and RFI Mitigation Implementations

Scope: which existing library has the best (fastest, most correct, most portable-licensed) implementation of each I/O and RFI function, as input to a unified Rust (CPU) / C++-CUDA (GPU) / Python stack. Research date: late Nov 2026. Status labels: items under "Cited Findings" were checked against a source in this session. Items under "Inferences" that come from background knowledge and were NOT re-checked are marked [unverified].

## Q1. Format readers/writers: which libraries handle filterbank, PSRFITS, PSRDADA, GUPPI/VDIF/baseband, and conversion between them, and how fast or robust are they?

### Takeaway
No single library covers every format well. `your` (GPL-3.0) is the broadest Python "unified reader" (read/write of fil and PSRFITS, plus PSRDADA writing). sigpyproc3 (MIT) is the best permissively licensed Python reader for fil and PSRFITS, with Numba kernels and 1/2/4-bit unpacking. PRESTO, PulsarX and TransientX read fil and PSRFITS in C/C++ but are GPL. DSPSR is still the de-facto reference for baseband formats (GUPPI raw, VDIF, DADA), but I did not re-verify its format list or license in this session.

### Cited Findings
- `your` ("Your Unified Reader") reads and writes SIGPROC filterbank and PSRFITS. For PSRDADA it can only write; reading is marked unsupported in its README. Licence: GPL-3.0. About 767 commits. — [thepetabyteproject/your GitHub](https://github.com/thepetabyteproject/your)
- The JOSS paper for `your` is Aggarwal et al. 2020, JOSS 5(55), 2750. The stated goal is to "unify the data processing across multiple commonly used formats". Latest version on PyPI is 0.6.x (0.6.7 per the search index). — [JOSS paper](https://joss.theoj.org/papers/10.21105/joss.02750.pdf); [PyPI your](https://pypi.org/project/your)
- `your` docs have module pages for Filwriter, Fitswriter, Psrdada, Psrfits, Pysigproc, plus utils pages "Rfi" and "Gpu". The index page itself lists no function-level detail. — [your readthedocs](https://your.readthedocs.io/en/latest/)
- sigpyproc3 reads SIGPROC filterbank, PSRFITS, time series and Fourier series. It uses Numba for its numerical kernels, which are JIT-compiled on the first call. Licence: MIT. 2.0.0 is the stable API; the last 0.x release was 0.5.5. — [FRBs/sigpyproc3 GitHub](https://github.com/FRBs/sigpyproc3)
- `sigpyproc.io.bits` unpacks 1-, 2- and 4-bit data stored in uint8 arrays, with a configurable bit order. Unpacking runs only for nbits in {1, 2, 4}; 8-bit data is used as is. — [sigpyproc3 io.bits docs](https://sigpyproc3.readthedocs.io/en/latest/generated/sigpyproc.io.bits.html)
- PRESTO v2's `readfile` and its general readers take SIGPROC filterbank and PSRFITS input. rfifind writes a mask that the downstream tools apply by replacing flagged samples with median values. — [PRESTO search tutorial (Ransom)](https://www.cv.nrao.edu/~sransom/PRESTO_search_tutorial.pdf)
- TransientX reads filterbank by default and PSRFITS via `--psrfits`. It depends on Boost >= 1.56, PlotX, XLibs, SOFA/ERFA and YMW16. Licence: GPL-3.0. — [ypmen/TransientX GitHub](https://github.com/ypmen/TransientX)
- PDAT (PulsarDataToolbox) supports only PSRFITS. — [PyPI pdat](https://pypi.python.org/pypi/pdat)
- FLITS is a browser-based burst-inspection tool, not a reader library. It reads SIGPROC .fil, PSRFITS search/fold (.fits/.sf) and CHIME/FRB HDF5. — [PyPI flits](https://pypi.org/project/flits/0.7.3/)

### Inferences
- [unverified] Common background knowledge, worth checking before relying on it:
  - DSPSR reads PSRDADA, GUPPI/PUPPI raw, VDIF, Mark5, LOFAR and many other baseband formats, and `digifil`/`digifits` convert baseband to fil/PSRFITS search mode.
  - blimpy (Breakthrough Listen, BSD-3) reads and writes .fil and BL HDF5 (.h5) and converts between them (`fil2h5`/`h52fil`).
  - PSRDADA itself is the C ring-buffer library (shared-memory HDU, 4096-byte ASCII header).
  - PulsarX shares TransientX's XLibs I/O layer and is also GPL-3.0.
- For a permissive Rust/C++ core, the safest starting points are:
  - MIT sigpyproc3 logic for fil/PSRFITS header parsing and bit unpacking.
  - A clean-room implementation of the DADA header format (simple ASCII key-value pairs).
  - The published PSRFITS definition (ATNF) for the binary-table layout.
  - Porting GPL readers (your, PRESTO, TransientX/PulsarX) would make the derived work GPL.
- Memory mapping and streaming:
  - sigpyproc3 documents block-wise reading (`read_block`/`read_plan`), but memmap was not stated on the pages fetched.
  - A Rust implementation using `memmap2` plus a block iterator would likely match or beat the Python readers. This is an engineering judgement, not a benchmark.
- 16- and 32-bit (float) filterbank need no unpacking, only dtype casting. Sub-byte unpacking (1/2/4-bit) and the PSRFITS scale/offset/weights application (DAT_SCL, DAT_OFFS, DAT_WTS) are where correctness bugs usually appear. A unified codebase should test against PRESTO/DSPSR outputs as the reference.

### Gaps
- No published head-to-head throughput benchmarks (GB/s) for fil/PSRFITS readers across your, sigpyproc3, PRESTO, PulsarX and blimpy were found.
- I did not verify these licences or format lists in this session:
  - DSPSR, PSRCHIVE, PRESTO and blimpy licences (believed: PRESTO GPL-2.0, blimpy BSD-3, DSPSR/PSRCHIVE GPL/AFL mix).
  - DSPSR's exact format list.
  - Whether `your`'s Psrdada module reads from live ring buffers.
- GUPPI raw and VDIF Python readers (e.g. baseband/`baseband-tasks`, `rawspec`) were not checked in this session.

## Q2. Are there Rust crates for FITS, filterbank, PSRFITS and DADA, and how mature are they?

### Takeaway
Rust has a usable SIGPROC filterbank crate (`sigproc_filterbank`, MIT/Apache-2.0, pure safe Rust, u1/u2/u4/u8/... via types). It also has three FITS options:
- `fitsio`: cfitsio FFI, the most mature and complete.
- `fitsio-pure`: pure Rust, reads and writes BINTABLE, Apache-2.0, up to about 3-10x slower than cfitsio in some cases.
- `fitsrs`: pure Rust, but binary tables are read-only.

No PSRFITS-specific or PSRDADA crate was found, so both need to be built on top of a FITS crate or from scratch.

### Cited Findings
- `sigproc_filterbank`:
  - Pure-Rust parser and serializer for SIGPROC filterbank files. It denies unsafe code and is dual-licensed MIT/Apache-2.0.
  - Bit depth comes from the pushed element type. Sub-byte data uses the `ux` crate (u1, u2, u4).
  - Latest version found is 0.4.0 (Apr 2024). Repo: github.com/kiranshila/sigproc_filterbank.
  - Sources: [docs.rs sigproc_filterbank](https://docs.rs/sigproc_filterbank); [lib.rs](https://www.lib.rs/crates/sigproc_filterbank); [docs.rs 0.3.0](https://docs.rs/crate/sigproc_filterbank/0.3.0)
- `fitsio` (simonrw/rust-fitsio): an FFI wrapper around cfitsio via `fitsio-sys`. It needs the C library and a C toolchain. — [rust-fitsio GitHub](https://github.com/simonrw/rust-fitsio); [fitsio-sys docs](https://docs.rs/crate/fitsio-sys/0.5.7/source/src/lib.rs)
- `fitsio-pure`, from its own README (self-reported; version 0.25.1 on the fetched page, while the search index showed 0.12.0):
  - Licence and portability: Apache-2.0, no C dependencies, no unsafe code, works on wasm32 and no_std.
  - Binary tables: read and write.
  - Tile compression: RICE_1 and GZIP_1 read/write; HCOMPRESS_1 and PLIO_1 read only.
  - Speed vs cfitsio: "within 3x on large images". Column writes at about 1M rows are around 10x slower.
  - Comparison table: `fitsrs` reads binary tables but cannot write them, has read-only headers, and only limited GZIP/RICE compression.
  - Source: [fitsio-pure README](https://docs.rs/crate/fitsio-pure/latest/source/README.md)
- Other pure-Rust FITS crates exist: `fits-io` (async Tokio I/O) and `rsfitsio` (a Rust rewrite of cfitsio). — [fits-io docs](https://docs.rs/crate/fits-io/0.2.0); [rsfitsio docs](https://docs.rs/crate/rsfitsio/latest)
- The crate named `sigproc` (0.1.0) is generic signal processing, not the SIGPROC file format. — [docs.rs sigproc](https://docs.rs/crate/sigproc/0.1.0/features)

### Inferences
- PSRFITS search mode is a FITS BINTABLE (SUBINT) holding a large byte-array DATA column plus DAT_SCL/DAT_OFFS/DAT_WTS/DAT_FREQ columns.
  - Needing binary-table write rules out `fitsrs` as the writer.
  - `fitsio` (cfitsio) is the most battle-tested choice.
  - `fitsio-pure` is attractive for a dependency-free build, but its roughly 10x column-write gap suggests writing the SUBINT data through a direct, row-aligned raw-byte path.
- For reading at full speed, the most robust approach is probably to parse the FITS headers with a crate, then memory-map the SUBINT data region and decode rows directly. FITS binary tables are fixed-width big-endian rows, so this is straightforward.
- `sigproc_filterbank` 0.4.0 has seen no releases since April 2024 (by the dates found). Treat it as a reference or a fork base rather than a hard dependency for multi-GB streaming. Whether it supports memory mapping or streaming was not verified.
- PSRDADA: no crate found. A Rust client would need FFI to libpsrdada (shared memory plus semaphores) or a reimplementation of the ipcbuf protocol. Reading `.dada` files only needs a 4096-byte ASCII header parser, which is easy.

### Gaps
- No crate was found for PSRDADA, GUPPI raw or VDIF. Searching crates.io again for "vdif", "dada", "guppi" and "psrfits" is recommended (not done in this session because of the tool-call budget).
- No download or maintenance metrics were collected for any crate.

## Q3. Which RFI algorithms does the literature best support for real-time search, which have GPU implementations, and what are the licences?

### Takeaway
For real-time search, the strongest-supported set is:
- IQRM for narrow-band channel masking: non-parametric, trend-robust, MeerKAT/Jodrell results, permissive Python reference implementation.
- Spectral kurtosis (SK) for non-Gaussian and transient-in-frequency RFI. It is cheap on GPU: CHIME runs it at about 0.7% of X-engine compute, and CHORD builds GPU SK masks every millisecond.
- Zero-DM filtering for broadband impulsive RFI. It is trivially parallel but suppresses low-DM signals.

PRESTO rfifind is still the offline benchmark, and recent papers report beating it. clfd (MIT) covers folded-data cleaning.

GPU codes include AstroAccelerate (GPL-3.0 per an aggregator), heimdall (built-in narrowband and 0-DM excision, licence unverified), the ARTS/Apertif RFIm library, CHIME kotekan SK, SPOTLIGHT VOLT/STRIPE, and ASTROFLOW.

### Cited Findings
**IQRM**
- Morello, Rajwade & Stappers, MNRAS 510 (2022), arXiv:2108.12434. It is a non-parametric outlier detector that is robust to trends in sequential data and is designed for real-time use with massively multibeam telescopes. — [arXiv 2108.12434](https://arxiv.org/abs/2108.12434v2); [ar5iv](https://ar5iv.labs.arxiv.org/html/2108.12434)
- Method, as given in the authors' RFI2022 slides:
  - Compute a per-channel statistic, such as standard deviation, over short blocks.
  - Take lagged differences to remove trends.
  - Apply Tukey's rule using a robust sigma derived from the IQR.
  - On MeerKAT and Jodrell Bank data, this gave roughly 3x more genuine pulses and about 3x fewer spurious candidates.
  - Source: [RFI2022 Morello slides](https://events.ecmwf.int/event/258/contributions/2894/attachments/1577/2816/RFI2022_Morello.pdf)
- IQRM is mainly a bright narrow-band masking algorithm and is not effective against long-term broadband RFI. — [FRB search review, arXiv 2405.20716](https://arxiv.org/pdf/2405.20716)
- Software: IQRM is registered in the ASCL (2310.009). — [ASCL 2310.009](https://ascl.net/2310.009)
- sigpyproc3 includes `sigpyproc.core.rfi.iqrm_mask` (defaults: threshold 3 sigma, radius 5), and an `RFIMask` class whose `apply_method` takes "iqrm" or "mad". — [sigpyproc3 rfi docs](https://sigpyproc3.readthedocs.io/en/latest/generated/sigpyproc.core.rfi.html); [sigpyproc3 core API](https://sigpyproc3.readthedocs.io/en/stable/api/core.html)

**Zero-DM**
- The zero-DM filter subtracts the DM=0 time series from every channel. This removes broadband RFI but can remove real low-DM signals. — [arXiv 2405.20716](https://arxiv.org/pdf/2405.20716)
- The original zero-DM paper is Eatough, Keane & Lyne 2009, "An Interference Removal Technique for Radio Pulsar Searches". — [arXiv 0901.3993](https://arxiv.org/pdf/0901.3993)

**Comparative evaluation**
- Hombal, Levin, Stappers et al. (arXiv:2601.08351, Jan 2026, accepted to RASTI) compare IQRM, SKF and ZDMF.
  - They inject pulses into simulated RFI environments and run a single-pulse search pipeline across brightness, width and DM.
  - Which method was best is not stated in the abstract.
  - Source: [arXiv 2601.08351](https://www.arxiv.org/abs/2601.08351)

**Spectral kurtosis on GPU**
- CHIME computes SK estimates on GPUs with OpenCL kernels in kotekan, at about 0.7% of X-engine compute. — [Taylor et al., SK for CHIME, arXiv 1808.10365](https://arxiv.org/pdf/1808.10365)
- CHORD's FRB backend computes per-channel variance and SK on GPU every ~1 ms to build a real-time mask before dedispersion. Details are deferred to a future paper. — [CHORD overview, arXiv 2607.09374](https://arxiv.org/pdf/2607.09374)
- A CSIR paper describes a PyCUDA polyphase filterbank plus SK detector on a Tesla T4 at 8.88 GSamples/s. This is spectrometer-side, not an FRB search. — [CSIR researchspace](https://researchspace.csir.co.za/bitstreams/302ef908-b116-45c8-be79-915eddd09ee5/download)

**Other GPU RFI codes**
- RFIm (Sclocco & Vohl, arXiv:2001.03389): an open-source library of real-time RFI mitigation routines for many-core accelerators, tunable for performance portability. It runs in the ARTS/Apertif pipeline (73 Gb/s) with "minimal impact on the total execution time" and considerably fewer false positives. The abstract gives neither the licence nor the repository URL. — [arXiv 2001.03389](https://arxiv.org/abs/2001.03389)
- Heimdall (GPU) has these CLI options: baseline smoothing length (default 2 s), RFI excision threshold (default 5), and flags to disable narrowband or 0-DM RFI excision. — [heimdall wiki Use](https://sourceforge.net/p/heimdall-astro/wiki/Use/); [heimdall project](https://sourceforge.net/p/heimdall-astro)
- Heimdall++ (2025) describes heimdall's RFI stage as multi-level (broadband plus narrowband), with separate baseline-removal and normalization stages. — [Heimdall++ arXiv 2512.00398](https://arxiv.org/pdf/2512.00398)
- AstroAccelerate:
  - A GPU (NVIDIA) package for real-time time-domain radio data.
  - Licence listed as GPL-3.0 by the ecosyste.ms aggregator, not by the repo's own LICENSE, so verify it.
  - Last Zenodo release found: v1.8.1, Nov 2020.
  - Sources: [Oxford AstroAccelerate](https://oerc.ox.ac.uk/research-groups/astroaccelerate/about); [ecosyste.ms](https://awesome.ecosyste.ms/projects/github.com%2FAstroAccelerateOrg%2Fastro-accelerate); [Zenodo](https://zenodo.org/records/4282748)
- SPOTLIGHT (GMRT; Wani, Roy, ..., Adamek, Armour et al., arXiv:2608.19102, Aug 2026):
  - Two stages: VOLT (antenna-level voltage filtering before beamforming, for broadband impulsive RFI) and STRIPE (statistical filtering of beamformed data, for narrowband RFI).
  - Benchmarked against rfifind: 98% fewer false detections and 2.7x S/N improvement on recovered pulses.
  - Source: [arXiv 2608.19102](https://arxiv.org/abs/2608.19102)
- ASTROFLOW (arXiv:2511.02328): a C++/CUDA single-pulse pipeline (RFI excision, dedispersion, YOLO detector) that runs more than 10x faster than real time on consumer GPUs. Its RFI method was not detailed in the snippet. — [arXiv 2511.02328](https://arxiv.org/abs/2511.02328v1)

**CPU / folded-data codes**
- TransientX RFI options:
  - `--zapthre` threshold zapping
  - Kadane filter (`-z kadaneF`) with `zdot`
  - `--clean`
  - `--iqr` IQR-based baseline/noise estimation, recommended when there are many strong bursts
  - Licence: GPL-3.0.
  - Source: [ypmen/TransientX](https://github.com/ypmen/TransientX)
- clfd ("clean folded data"):
  - Built on PSRCHIVE and developed for the HTRU reprocessing.
  - Main step: reduce each (channel, subint) profile to features (ptp, std, var, lfamp, skew, kurtosis, acf), flag outliers with Tukey's rule, and set PSRCHIVE weights to 0.
  - Optional `--despike` step: time-phase spike subtraction.
  - Requires PSRCHIVE Python bindings. Licence: MIT.
  - Source: [v-morello/clfd](https://github.com/v-morello/clfd)
- PRESTO rfifind:
  - Uses time-domain plus frequency-domain block statistics to flag strong narrow-band and/or short-duration broadband RFI, and writes a mask. — [PRESTO tutorial](https://www.cv.nrao.edu/~sransom/PRESTO_search_tutorial.pdf)
  - It is still "the standard RFI mitigation tool for major FAST pulsar surveys" and is used as a benchmark. — [arXiv 2512.16278](https://arxiv.org/pdf/2512.16278)
  - A multi-beam eigenvalue-decomposition method reports removing more RFI than rfifind. A 2026 transient-search paper reports +17.6% median folded S/N over rfifind and says rfifind sometimes corrupted input data. These are per the search summary; the claims were not individually verified. — [arXiv 2512.16278](https://arxiv.org/pdf/2512.16278); [SPOTLIGHT arXiv 2608.19102](https://arxiv.org/pdf/2608.19102)

### Inferences
- Recommended "best implementation to port" per function. These are judgements from the evidence above.

  | Function | Recommended source |
  |---|---|
  | IQRM | Morello's reference Python (permissive per ASCL entry; licence not re-verified) or sigpyproc3 `iqrm_mask` (MIT) |
  | Zero-DM | Trivial; write from the Eatough et al. 2009 description |
  | SK / generalized SK | Write from the Nita & Gary papers [unverified citation: Nita & Gary 2010 PASP 122, 595 and MNRAS 406, L60]. Use CHIME's kotekan kernels as a GPU design reference; their licence was not checked |
  | Folded cleaning | clfd (MIT) feature + Tukey approach |
  | rfifind-style block statistics | Reimplement from the description; PRESTO is GPL [unverified] |
  | Kadane / zdot (TransientX) | GPL-3.0, so a clean-room reimplementation is needed for a permissive codebase |

- Real-time GPU ordering seen in deployed systems (CHORD, heimdall, SPOTLIGHT): per-block channel statistics (variance/SK/IQRM) → channel mask → zero-DM or broadband time-sample excision → baseline removal/normalization → dedispersion.
  - All of these are embarrassingly parallel reductions over (channel, block), well suited to CUDA kernels. A Rust CPU version can use SIMD plus rayon.
- `your` exposes an RFI utils module (savgol / SK-based channel flagging per [unverified] background knowledge) and a GPU module (likely CuPy). Being GPL-3.0, it should be used as a validation oracle rather than ported.
- [unverified] Further background points, worth checking before relying on them:
  - PSRCHIVE `paz` (zapping by channel/subint/median filter) and `pazi` (interactive) work only on folded archives.
  - DSPSR has built-in SK (`-skz`, generalized SK per Nita & Gary) on baseband/filterbank data.
  - The MeerTRAP pipeline (Mosaic/AstroAccelerate-based) uses IQRM-like and zero-DM cleaning.
  - Breakthrough Listen's rawspec/hyperseti do not focus on pulsar RFI.
  - ASKAP CRAFT uses a CUDA FREDDA search with its own channel-kurtosis/zero-DM flagging.

### Gaps
- No quantitative head-to-head results were found for IQRM vs SK vs zero-DM; the full text of arXiv 2601.08351 is needed.
- Licences not confirmed in this session: heimdall, RFIm, kotekan, PRESTO, PSRCHIVE, DSPSR, blimpy, and the IQRM reference repo. AstroAccelerate GPL-3.0 comes only from an aggregator.
- Not covered by fetched sources: MeerTRAP, ASKAP CRAFT (FREDDA) and Breakthrough Listen GPU RFI specifics; DSPSR SK implementation details; the generalized SK (Nita & Gary) primary sources.
