# Math maps

One file per source library. Each map lists, in pipeline order, every mathematical operation that library performs: the equation with its symbols and units, where it lives in the library's source, the paper it comes from, and which of our modules does the same thing (or "planned"). Equations not checked against a primary source are marked **[unverified]**. The same equations are gathered by topic in [05-algorithms-math.md](../reference/05-algorithms-math.md); what our code has already taken is recorded in [CREDITS.md](../CREDITS.md); the reason these maps exist is set out in [06-our-plan.md](../reference/06-our-plan.md).

| Library | Role | Map |
|---|---|---|
| PRESTO | CPU search toolkit: rfifind, DDplan, prepsubband, accelsearch, coherentsearch, prepfold | [presto.md](presto.md) |
| PrestoZL | GPU fork of PRESTO's accelsearch (jerk) and prepsubband | [prestozl.md](prestozl.md) |
| dmt | Dedispersion library: FDMT, direct, subband, coherent, with exact per-trial noise | [dmt.md](dmt.md) |
| AstroAccelerate | GPU library: direct dedispersion, SPDT single-pulse search, FDAS, jerk | [astroaccelerate.md](astroaccelerate.md) |
| heimdall | Real-time GPU single-pulse search: baseline, prefix-sum boxcars, friends-of-friends | [heimdall.md](heimdall.md) |
| TransientX | Real-time CPU single-pulse search: RFI chain, subband dedispersion, DBSCAN | [transientx.md](transientx.md) |
| PulsarX | Fold-and-optimise back end: pruned FDMT, Tikhonov folding, DM/F0/F1 grid | [pulsarx.md](pulsarx.md) |
| riptide | FFA periodicity search with exact boxcar S/N | [riptide.md](riptide.md) |
| peasoup | GPU binary-pulsar FFT search with time-domain acceleration resampling | [peasoup.md](peasoup.md) |
| sigpyproc3 | Python filterbank toolbox: readers, ZERO_OFF rules, IQRM, reference searches | [sigpyproc3.md](sigpyproc3.md) |
| SIGPROC | Filterbank and time-series format: header encoding, ID tables, LSB-first packing | [sigproc.md](sigproc.md) |
| PSRDADA | Shared-memory ring buffers and `.dada` header for real-time data | [psrdada.md](psrdada.md) |
| DSPSR | Baseband processing: coherent dedispersion, spectral kurtosis, detection, folding | [dspsr.md](dspsr.md) |
| PSRCHIVE | Folded-archive tools and `pat` PGS template-matching TOAs | [psrchive.md](psrchive.md) |
| your | Unified reader, format writer and FETCH cut-out maker | [your.md](your.md) |
| FETCH | Two-branch CNN classifier for single-pulse candidates | [fetch.md](fetch.md) |
| TEMPO2 | Barycentring, timing-model fits, polyco and Chebyshev predictors | [tempo2.md](tempo2.md) |
| PINT | Python timing package: delay chain, polycos, WLS and GLS fitters | [pint.md](pint.md) |
| clfd | Folded-data RFI cleaning by profile features and Tukey fences | [clfd.md](clfd.md) |

## Conventions used in every map

- Math is GitHub math: inline `$...$` and fenced `math` blocks. Frequencies are in MHz and times in seconds unless a map says otherwise.
- Tags: [V] read in the primary paper or specification; [C] transcribed from the library's source by the reference set; [D] derived here from tagged equations; **[unverified]** not checked against a primary source.
- "Ours" links point into `filterbank/` with relative paths; "planned" means the module is in [06-our-plan.md](../reference/06-our-plan.md) but not yet written.
- Each map has one Mermaid flowchart showing the library's processing order.
