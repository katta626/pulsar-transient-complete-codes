# On-disk data formats for radio pulsar and transient searching (SIGPROC, PSRFITS, PSRDADA, PRESTO, heimdall/FETCH)

Method note: the details below come mostly from reading the source code (shallow clones made 2026-11-24) and the ATNF PSRFITS v6.1 documentation page. Source snapshots used:

| Repo | URL | Commit (date) |
|---|---|---|
| SIGPROC (SixByNine fork, the maintained one) | https://github.com/SixByNine/sigproc | c8bcbd7 (2026-05-26) |
| PRESTO | https://github.com/scottransom/presto | 4a5bd80 (2026-11-18) |
| sigpyproc3 | https://github.com/FRBs/sigpyproc3 | 5d8d38f (2026-10-17) |
| your | https://github.com/thepetabyteproject/your | 5126d3b (2026-09-19) |
| heimdall | https://sourceforge.net/p/heimdall-astro/code/ | 69a33da (2023-06-15) |
| PSRDADA | https://sourceforge.net/p/psrdada/code/ | 151198e (2026-04-02) |
| PSRFITS definition v6.1 | https://www.atnf.csiro.au/research/pulsar/psrfits_definition/PsrfitsDocumentation.html (index: https://www.atnf.csiro.au/research/pulsar/psrfits_definition/Psrfits.html) | page fetched 2026-11-24 |

Line numbers like `read_header.c:L..` refer to these snapshots. GitHub file URLs use the form `https://github.com/<repo>/blob/master/<path>` (sigpyproc3/your use `main`).

---

## 1. SIGPROC filterbank (.fil) and time series (.tim)

### Takeaway
A SIGPROC file is a binary keyword header (each keyword is `int32 length` + ASCII name, then a typed value in the writer's native byte order) bracketed by `HEADER_START`/`HEADER_END`, followed by raw samples in time-major order (one full spectrum per time sample, channel 1 = `fch1` first). Sub-byte samples are packed least-significant-bits-first in SIGPROC itself (the opposite of PSRFITS). `.tim` is the same header with `data_type=2`, `nchans=1`, followed by float32 samples.

### Cited Findings

**Keyword encoding**
- Strings are written by `send_string()`: `int len = strlen(string); fwrite(&len, sizeof(int),1); fwrite(string, 1, len)`. No NUL terminator is written. Typed values are written by `send_int` (C `int`, 4 B), `send_double` (8 B), `send_float` (4 B), `send_char` (1 B) and `send_long` (C `long`, 8 B on LP64), each as `send_string(name)` followed by the raw value. Byte swapping happens only when the global `swapout` is set — [sigproc src/send_stuff.c](https://github.com/SixByNine/sigproc/blob/master/src/send_stuff.c)
- String-valued keywords (`rawdatafile`, `source_name`) are written as two consecutive strings: the keyword, then the value string. The reader sets an `expecting_rawdatafile`/`expecting_source_name` flag and takes the next string as the value — [sigproc src/read_header.c](https://github.com/SixByNine/sigproc/blob/master/src/read_header.c)
- `get_string()` reads an `int nchar` and gives up (returns "ERROR") if `nchar>80 || nchar<1`. The header must start with `HEADER_START`, otherwise `read_header` rewinds and returns 0 ("not in standard format"). An unknown keyword is a fatal error (`read_header - unknown parameter`). The return value is the total header length in bytes, which is also checked against `ftell` — [sigproc src/read_header.c](https://github.com/SixByNine/sigproc/blob/master/src/read_header.c)
- PRESTO's `read_filterbank_header()` follows the same algorithm: it also exits on an unknown keyword, sets `ibeam=1` by default, and computes `N = (filelen - headerlen) / (nchans / (8/nbits))` for nbits ≤ 8 — [presto src/sigproc_fb.c](https://github.com/scottransom/presto/blob/master/src/sigproc_fb.c)

**Keyword table (types as read by `read_header.c`, globals declared in `header.h`)**

| Keyword | C type on disk | Meaning / units | Source notes |
|---|---|---|---|
| `HEADER_START` / `HEADER_END` | string only | header delimiters | read_header.c |
| `rawdatafile` | string (next string) | original raw file name, ≤80 chars | `char rawdatafile[80]` |
| `source_name` | string (next string) | source name, ≤80 chars | `char source_name[80]` |
| `telescope_id` | int32 | telescope code (table below) | aliases.c |
| `machine_id` | int32 | backend code (table below) | aliases.c |
| `data_type` | int32 | 0 raw data, 1 filterbank, 2 time series, 3 pulse profiles, 4 amplitude spectrum, 5 complex spectrum, 6 dedispersed subbands | `data_category()` in aliases.c; sigpyproc adds 10 = "PSRFITS" |
| `barycentric` | int32 | 1 if times are barycentric | added "Aug 20, 2004 DRL" per header.h comment |
| `pulsarcentric` | int32 | 1 if times are pulsar-centric | same |
| `az_start` | double | telescope azimuth at start, degrees | |
| `za_start` | double | zenith angle at start, degrees | |
| `src_raj` | double | J2000 RA packed as **hhmmss.s** (e.g. 053431.97 = 05:34:31.97) | PRESTO decodes `h=floor(raj/10000)`, `m=floor((raj-h*10000)/100)`, `s=raj-h*10000-m*100` and rejects values outside [0, 240000) |
| `src_dej` | double | J2000 Dec packed as **±ddmmss.s** | PRESTO accepts \|dej\| ≤ 900000; the sign belongs to the whole value, so −00:MM:SS needs care (comment in sigproc_fb.c) |
| `tstart` | double | MJD (UTC, topocentric unless `barycentric`) of the first sample | |
| `tsamp` | double | sampling interval, **seconds** | |
| `nbits` | int32 | bits per sample (1, 2, 4, 8, 16, 32) | |
| `nsamples` | int32 | number of samples; read **only for backwards compatibility** and discarded (`itmp`) | read_header.c comment |
| `fch1` | double | centre frequency of the first channel in the file, MHz | |
| `foff` | double | channel offset, MHz; negative means frequency decreases with channel index | |
| `fchannel` | double (repeated) | explicit per-channel frequency table between `FREQUENCY_START`/`FREQUENCY_END`; when used, `fch1=foff=0.0` | `double frequency_table[4096]` (hard limit of 4096 channels) |
| `nchans` | int32 | number of channels | |
| `nifs` | int32 | number of IFs/polarisations | |
| `refdm` | double | reference DM (pc cm⁻³), used for `.tim`/subbands | |
| `period` | double | folding period, seconds (profiles, data_type 3) | |
| `nbins` | int32 | number of phase bins (data_type 3) | |
| `npuls` | C `long` (8 B on LP64) | number of pulses (binary profile format) | read as `sizeof(npuls)` |
| `nbeams` | int32 | number of beams | |
| `ibeam` | int32 | beam number | |
| `signed` | char (1 B) | signedness of 8-bit data; compared against compile-time `OSIGN` | `#define OSIGN 1` (unsigned build) |

Sources: [sigproc src/header.h](https://github.com/SixByNine/sigproc/blob/master/src/header.h), [sigproc src/read_header.c](https://github.com/SixByNine/sigproc/blob/master/src/read_header.c), [sigproc src/aliases.c](https://github.com/SixByNine/sigproc/blob/master/src/aliases.c), [presto src/sigproc_fb.c](https://github.com/scottransom/presto/blob/master/src/sigproc_fb.c)

- sigpyproc3 keeps a larger key table (`header_keys` in `io/params.py`). On top of the SIGPROC keys it reads `flux` (d), `pb` (d), `ecc` (d), `asini` (d), `accel` (d), `bandwidth`/`fbottom`/`ftop` (d), `obs_date`/`obs_time` (str), `hdrlen`/`sampsize`/`orig_hdrlen`/`new_hdrlen` (uint) and `signed` ("b"). Its integers are struct code "I" (uint32). Units given there: az/za in degrees, src_raj "(HH:MM:SS.ss)", tstart MJD, tsamp s, fch1/foff MHz, refdm pc cm⁻³, period s, pb hours — [sigpyproc3 src/sigpyproc/io/params.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/io/params.py)

**Telescope and machine IDs (`aliases.c`)**
- telescope_id: 0 Fake, 1 Arecibo, 2 Ooty, 3 Nancay, 4 Parkes, 5 Jodrell, 6 GBT, 7 GMRT, 8 Effelsberg, 9 140ft, 10 SRT, 64 MeerKAT, 65 KAT-7, 82 eMerlin. TEMPO site codes: 1→'3', 3→'f', 4→'7', 5→'8', 6→'1', 8→'g', 9→'a', 10→'z', 64→'m', 65→'k' — [sigproc aliases.c](https://github.com/SixByNine/sigproc/blob/master/src/aliases.c)
- machine_id: 0 FAKE, 1 PSPM, 2 WAPP, 3 AOFTM, 4 BPP, 5 OOTY, 6 SCAMP, 7 GMRTFB, 8 PULSAR2000, 9 PARSPEC, 10 BPSR, 14 GMRTNEW, 64 KAT, 65 KAT-DC2, 82 loft-e, 83 ROACH, 88 AT-PAF — [sigproc aliases.c](https://github.com/SixByNine/sigproc/blob/master/src/aliases.c)
- PRESTO keeps its own, partly different ID tables in `get_telescope_name()`/`get_backend_name()` (cases include 11, 12, 20, 21, 30, 31, 32) — [presto src/sigproc_fb.c](https://github.com/scottransom/presto/blob/master/src/sigproc_fb.c). The two tables do not fully match, so IDs above 10 depend on which tool wrote the file.

**Writer behaviour (`filterbank_header.c`)**
- Write order: `HEADER_START`, `rawdatafile`+name, optional `source_name`, `machine_id`, `telescope_id`, coordinates via `send_coords`, then for a time series `data_type=2`, `refdm` and `nchans=1`, or for a filterbank `data_type=1`, `fch1`, `foff`, `nchans`. Then `nbeams`, `ibeam`, `nbits`, `tstart` (+ start_time/86400), `tsamp`, `nifs` and `HEADER_END`. When the band is inverted on output it writes `fch1 = fch1 + foff*(nchans-1)` and `foff = -foff` — [sigproc src/filterbank_header.c](https://github.com/SixByNine/sigproc/blob/master/src/filterbank_header.c)
- Bug in `send_coords()`: the guard `if ((raj != 0.0) || (raj != -1.0))` is always true, so `src_raj`, `src_dej`, `az_start` and `za_start` are always written, even when they hold the sentinels 0 or −1 — [sigproc send_stuff.c](https://github.com/SixByNine/sigproc/blob/master/src/send_stuff.c)

**Data ordering and sub-byte packing**
- `read_block()` handles nbits 1/2/4/8/16/32 by unpacking a flat byte stream into floats in file order. For 1-bit: `block[k++]=charblock[i]&1; charblock[i]>>=1;` (LSB first). For 2-bit, `char2fourints` takes i = bits 0–1, j = bits 2–3, k = bits 4–5, l = bits 6–7. For 4-bit, `char2ints` takes the low nibble first (the comment says "first integer passed (i) takes the lower four bits"). So **SIGPROC packs the earliest channel in the least-significant bits** — [sigproc src/pack_unpack.c](https://github.com/SixByNine/sigproc/blob/master/src/pack_unpack.c), [read_block.c](https://github.com/SixByNine/sigproc/blob/master/src/read_block.c)
- sigpyproc3 defaults to `bitorder="little"` for 1/2/4-bit filterbank (`default_bitorder = {1:"little",2:"little",4:"little"}`), but uses `bitorder="big"` in `_expand_subbyte()` for PSRFITS ("shared big-endian unpacker") — [sigpyproc3 io/bits.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/io/bits.py), [io/pfits.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/io/pfits.py)
- PRESTO `convert_filterbank_block()` **always reverses channel order** while reading (for 32/16/8-bit, `outdata[ii] = data[nchan-1-ii]`). For packed data it walks the bytes from the end of the spectrum and emits MSB-first within each byte. Combined, these give a correct low-to-high frequency ordering of LSB-first-packed data. It sets `df=fabs(foff)`, `lo_freq = fch1 - (nchan-1)*df`, `hi_freq=fch1`, and `apply_flipband=0` "automatically" only when `foff<0` — [presto src/sigproc_fb.c](https://github.com/scottransom/presto/blob/master/src/sigproc_fb.c)
- Header-derived values PRESTO uses from a .fil: nbits, signed, nchans, nifs (if `nifs>1` it sets `sumifs=0` and `num_polns=nifs`), tsamp, foff, fch1, az_start, za_start, ibeam, telescope_id, machine_id, tstart, src_raj and src_dej. It ignores `nsamples` and computes N from the file size — [presto src/sigproc_fb.c](https://github.com/scottransom/presto/blob/master/src/sigproc_fb.c)
- heimdall's `SigprocFile` uses only `nchans`, `nbits`, `ibeam`, `tsamp` (seconds), `fch1` (f0), `foff` (df) and `tstart` (converted to `utc_start`) — [heimdall Formats/SigprocFile.C](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Formats/SigprocFile.C)
- `your`'s `pysigproc` maps only nbits 8→uint8, 16→uint16 and 32→float32, and raises `nbits=%d not supported` for anything else. It does not handle sub-byte filterbanks — [your formats/pysigproc.py](https://github.com/thepetabyteproject/your/blob/main/your/formats/pysigproc.py)

**.tim**
- For a time series the writer emits `data_type=2`, `refdm`, `nchans=1` and `fch1` (reference frequency); `header.c` prints "Reference DM" and "Reference frequency" for data_type 2 — [sigproc filterbank_header.c](https://github.com/SixByNine/sigproc/blob/master/src/filterbank_header.c), [header.c](https://github.com/SixByNine/sigproc/blob/master/src/header.c)

### Inferences
- The file has no endianness marker. Values are written in the writer's native order unless `swapout` is set. In practice all modern files are little-endian, and readers that see an implausible `int32` string length can try byte-swapping.
- PRESTO's `.fil` reader hardcodes the assumption `foff<0` (fch1 = highest frequency). For `foff>0` files the always-reverse path plus `lo_freq = fch1-(n-1)|foff|` gives a wrong frequency axis unless flipping is handled elsewhere. Writers should prefer `foff<0`, or verify with PRESTO's `readfile`.
- A 1/2/4-bit file can unpack into a byte-wise channel permutation if the packing convention is misread (SIGPROC LSB-first vs PSRFITS MSB-first). Because the permutation repeats every 2/4/8 channels, it is easy to miss on a plot.

### Gaps
- I did not find an original 1990s "format spec" document beyond the sigproc source and its LaTeX manual (`src/documentation.tex`, not read in full). The sigproc.sourceforge.net pages were not fetched.
- PRESTO's handling of `foff>0` was inferred from code reading and was not tested.

---

## 2. PSRFITS (ATNF/CSIRO definition, Hotan, van Straten & Manchester 2004, PASA 21, 302)

### Takeaway
PSRFITS is a FITS file with a primary HDU that describes the observation plus a set of binary-table HDUs. The `SUBINT` table holds the data. In search mode each row is `NSBLK` samples × `NPOL` × `NCHAN` as a byte array, sub-byte samples packed **earliest sample in the highest-order bits**. Physical values are `(DATA − ZERO_OFF)·DAT_SCL + DAT_OFFS`, optionally multiplied by `DAT_WTS`.

### Cited Findings

**Primary HDU keywords (v6.1)** — all from [ATNF PsrfitsDocumentation](https://www.atnf.csiro.au/research/pulsar/psrfits_definition/PsrfitsDocumentation.html)

| Keyword | Type/unit | Meaning |
|---|---|---|
| HDRVER, FITSTYPE | str | '6.1', 'PSRFITS' |
| DATE | str | file creation UTC (YYYY-MM-DDThh:mm:ss) |
| OBSERVER, PROJID, TELESCOP | str | TELESCOP "should conform with TEMPO2 usage" |
| ANT_X/Y/Z | D, m | antenna ITRF coordinates |
| FRONTEND, IBEAM | str | receiver ID, beam ID (beam 0 normally on-axis) |
| NRCVR | int | receiver polarisation channels (normally 2) |
| FD_POLN | str | 'LIN' or 'CIRC' |
| FD_HAND, FD_SANG, FD_XYPH | int, deg, deg | feed handedness (±1), feed angle of E-vector, cal cross-phase |
| BACKEND, BECONFIG | str | backend name, configuration file |
| BE_PHASE | int | 0/+1/−1 cross-product phase convention (sign flip conjugates the field, flips Stokes V) |
| BE_DCC | int | lower-sideband conjugation corrected flag |
| BE_DELAY | s | backend propagation delay (positive if time tag is early) |
| TCYCLE | s | on-line correlator cycle time |
| OBS_MODE | str | 'PSR' (fold), 'CAL' (fold on cal), 'SEARCH' (streamed multichannel) |
| DATE-OBS | str | approximate UTC start |
| OBSFREQ | MHz | "Radio frequency at the mid-point of channel OBSNCHAN/2 of the original observed spectrum, where channels are numbered from 1 to OBSNCHAN" |
| OBSBW | MHz | OBSNCHAN × channel bandwidth (nominal total) |
| OBSNCHAN | int | channels in the original spectrum |
| CHAN_DM | pc cm⁻³ | DM used for on-line (normally coherent) dedispersion |
| PNT_ID, SRC_NAME | str | pointing label, source name |
| COORD_MD, EQUINOX | str, float | J2000/GALACTIC/ECLIPTIC; equinox is a float |
| RA, DEC | str | 'hh:mm:ss.ssss', '−dd:mm:ss.sss' |
| BMAJ, BMIN, BPA | deg | beam |
| STT_CRD1/2, STP_CRD1/2, TRK_MODE, SCANLEN | | start/stop coordinates, TRACK/SCANGC/SCANLAT, requested length (s) |
| FD_MODE, FA_REQ | | feed tracking (FA/CPA/SPA/TPA), requested angle |
| CAL_MODE, CAL_FREQ, CAL_DCYC, CAL_PHS, CAL_NPHS | | injected-cal settings |
| STT_IMJD | int (J) | integer UTC MJD of start. Start = leading edge of first bin (fold) or first sample (search) |
| STT_SMJD | int (J), s | integer seconds past 0h UTC |
| STT_OFFS | D, s | fractional seconds |
| STT_LST | D, s | approximate LST at start |

**SUBINT HDU keywords** — [ATNF PsrfitsDocumentation](https://www.atnf.csiro.au/research/pulsar/psrfits_definition/PsrfitsDocumentation.html)

| Keyword | Meaning |
|---|---|
| EPOCHS | VALID / MIDTIME / STT_MJD (default STT_MJD if absent): how subint epochs relate to phase |
| INT_TYPE, INT_UNIT | time axis type (TIME, BINPHSPERI, BINLNGASC…), unit (SEC, PHS, DEG) |
| SCALE | intensity units (FluxDen/RefFlux/Jansky) |
| POL_TYPE | e.g. 'AA+BB' (NPOL=1), 'AABB' (NPOL=2), 'AABBCRCI', 'IQUV' |
| NPOL | number of polarisation products in DATA |
| TBIN | s; bin width (fold) or sampling interval (search) |
| NBIN | profile bins (fold), normally 1 for search |
| NBIN_PRD, PHS_OFFS | gated-data bins per period, phase of bin 0 |
| NBITS | bits per datum (search); ignored for fold |
| ZERO_OFF | "offset of the nominal data zero in digitiser units (actual data value = unsigned integer value − ZERO_OFF)" |
| SIGNINT | 1 for signed integers in search mode, else 0 |
| NSUBOFFS | accumulated subint count at start of this file for time-split files (true index = file index + NSUBOFFS) |
| NCHAN | channels in this file |
| CHAN_BW | MHz; "CHAN_BW < 0 implies decreasing frequency with increasing channel number" |
| DM, RM | DM/RM used for post-detection dedispersion / de-Faraday |
| NCHNOFFS | true channel index of first channel for frequency-split files |
| NSBLK | samples per row (search; else 1); "typically 4096" |
| NSTOT | total valid samples in this file (last row may be partly filled); per-file for split files |

**SUBINT columns**

| Column | TFORM | Unit | Meaning |
|---|---|---|---|
| INDEXVAL | 1D | | time-like coordinate if INT_TYPE ≠ TIME |
| TSUBINT | 1D | s | duration of subint/row |
| OFFS_SUB | 1D | s | time since observation start **at the centre** of each row |
| LST_SUB | 1D | s | LST at row centre |
| RA_SUB, DEC_SUB | 1D | deg | J2000 pointing at row centre |
| GLON_SUB, GLAT_SUB | 1D | deg | Galactic pointing |
| FD_ANG, POS_ANG, PAR_ANG | 1E | deg | feed, position and parallactic angles |
| TEL_AZ, TEL_ZEN | 1E | deg | azimuth, zenith angle at row centre |
| AUX_DM, AUX_RM | 1D | | additional time-varying DM/RM |
| DAT_FREQ | NCHAN D | MHz | centre frequency of each channel (sub-band = weighted mean) |
| DAT_WTS | NCHAN E | | weights 0–1 (0 = RFI or outside band) |
| DAT_OFFS | NCHAN·NPOL E | | "value subtracted from data to give zero mean for each channel" |
| DAT_SCL | NCHAN·NPOL E | | scale "(outval=dataval*scl + offs)" |
| DATA | I (fold) or B (search) | Jy | TDIM (NBIN,NCHAN,NPOL) fold or (NCHAN,NPOL,NSBLK·NBITS/8) search |

Source: [ATNF PsrfitsDocumentation](https://www.atnf.csiro.au/research/pulsar/psrfits_definition/PsrfitsDocumentation.html)

**Scaling equation, packing, fold vs search**
- Fold mode: "16-bit signed integers with elements of the data array in bin, channel and polarisation order with the pulse profile bins in contiguous locations", and "Real value = DATA value * DAT_SCL + DAT_OFFS" — [ATNF doc](https://www.atnf.csiro.au/research/pulsar/psrfits_definition/PsrfitsDocumentation.html)
- Search mode: 1/2/4/8-bit signed or unsigned integers written as a byte array, "packed with earlier samples in higher-order bits of the byte (i.e., 'big-endian')". Element order is channel, polarisation, sample, with channels contiguous — [ATNF doc](https://www.atnf.csiro.au/research/pulsar/psrfits_definition/PsrfitsDocumentation.html)
- Unsigned data: "Normally ZERO_OFFS = 2^(NBIT − 1) − 0.5 … Real value = (DATA value − ZERO_OFFS) * DAT_SCL + DAT_OFFS". The prose spells it ZERO_OFFS but the keyword is `ZERO_OFF` — [ATNF doc](https://www.atnf.csiro.au/research/pulsar/psrfits_definition/PsrfitsDocumentation.html)
- Search spectra are usually normalised as (S−R)/R before few-bit digitisation, and "are normally analysed directly without application of the scale factors and offsets". Time-split files are assumed to sample continuously across files — [ATNF doc](https://www.atnf.csiro.au/research/pulsar/psrfits_definition/PsrfitsDocumentation.html)

**Other HDUs (EXTNAMEs in v6.1)**
- HISTORY (processing history: DATE_PRO, PROC_CMD, SCALE, POL_TYPE, NSUB, NPOL, NBIN, NBIN_PRD, TBIN, CTR_FREQ, NCHAN, CHAN_BW, DM, RM, PR_CORR, FD_CORR, BE_CORR, RM_CORR, DEDISP, DDS_MTHD, SC_MTHD, CAL_MTHD, CAL_FILE, RFI_MTHD, RM_MODEL, AUX_RM_C, DM_MODEL, AUX_DM_C). Also OBSDESCR (free text), PSRPARAM (ephemeris text lines), POLYCO (DATE_PRO, POLYVER, NSPAN [min], NCOEF ≤15, NPBLK, NSITE, REF_FREQ, PRED_PHS, REF_MJD, REF_PHS, REF_F0, LGFITERR, COEFF), T2PREDICT (Tempo2 predictor text), COHDDISP (coherent dedispersion parameters), BANDPASS (original bandpasses; "Spectral value in Janskys = DATA value * DAT_SCL + DAT_OFFS"), FLUX_CAL, CAL_POLN, FEEDPAR, SPECKURT, SUBINT, DIG_STAT, DIG_CNTS — [ATNF doc](https://www.atnf.csiro.au/research/pulsar/psrfits_definition/PsrfitsDocumentation.html)

**How tools read PSRFITS search data**
- PRESTO `read_PSRFITS_files()` reads CHAN_DM, STT_IMJD/SMJD/OFFS, TBIN, NPOL, POL_TYPE, NSBLK, NSUBOFFS and ZERO_OFF (`s->zero_offset = fabs(...)`, so a negative ZERO_OFF is silently made positive). It also reads OFFS_SUB to detect missing rows ("NSUBOFFS reports %d previous rows but OFFS_SUB implies %d. Using OFFS_SUB"). If all OFFS_SUB are 0 it assumes no missing rows; if only some are 0 that is an error ("Not good PSRFITS") — [presto src/psrfits.c](https://github.com/scottransom/presto/blob/master/src/psrfits.c)
- PRESTO applies DAT_WTS, DAT_OFFS and DAT_SCL only when the first row is non-trivial (weights ≠1, offsets ≠0, scales ≠1), unless the user forces it. The formula is `((raw − zero_offset)*scale + offset)*weight`, applied to every nbits including 16- and 32-bit — [presto src/psrfits.c](https://github.com/scottransom/presto/blob/master/src/psrfits.c)
- PRESTO computes `df` from DAT_FREQ endpoints, `(f[N−1]−f[0])/(N−1)`, and only warns ("Channel spacing changes") if spacing varies by >1e-7 MHz. It does not support non-uniform channels. If `hi_freq < lo_freq` it sets `apply_flipband=1` and negates df — [presto src/psrfits.c](https://github.com/scottransom/presto/blob/master/src/psrfits.c)
- PRESTO unpacks PSRFITS MSB-first (`cdata[jj]=uctmp>>4; cdata[jj+1]=uctmp&0x0F` for 4-bit, and `>>6, >>4, >>2, &3` for 2-bit). It byte-reverses 8-byte groups for 1-bit data only for the special case TELESCOP=Parkes and BACKEND=FB_1BIT ("Flipping bit ordering since Parkes FB_1BIT data"). With NPOL>1 it sums polarisations if POL_TYPE starts "AABB" or NPOL==2, otherwise it uses one chosen pol — [presto src/psrfits.c](https://github.com/scottransom/presto/blob/master/src/psrfits.c)
- sigpyproc3 `PFITSFile.read_subint_pol()` follows "dspsr unpacker rules". ZERO_OFF is subtracted except for **2-bit** ("2-bit levels are already the dspsr lookup") and **16-bit** ("16-bit is a signed integer and dspsr does not subtract ZERO_OFF from it"). It returns Stokes I in DAT_FREQ order without reversing: PPQQ/Coherence give P0+P1, Intensity/Stokes give product 0. It refuses files whose OFFS_SUB spacing implies gaps ("split files and row gaps are not supported") and requires uniform DAT_FREQ (`_require_uniform_freqs`). `tstart = STT_IMJD + (STT_SMJD + STT_OFFS) s`, with NSUBOFFS not added. Its PSRFITS→SIGPROC key map is IBEAM→ibeam, NBITS→nbits, OBSNCHAN→nchans, SRC_NAME→source_name, RA→src_raj, DEC→src_dej, CHAN_BW→foff, TBIN→tsamp — [sigpyproc3 io/pfits.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/io/pfits.py), [io/params.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/io/params.py)
- `your` `PsrfitsFile.get_data()` applies `data *= DAT_SCL; data += DAT_OFFS; data *= DAT_WTS` and **never reads ZERO_OFF** (no occurrence in the package). Its 2/4-bit unpackers are MSB-first. For AABB/AABBCRCI it averages AA and BB (×0.5). For IQUV it can form I, RCP/LCP ((I±V)/2) or linear ((I±Q)/2) — [your formats/psrfits.py](https://github.com/thepetabyteproject/your/blob/main/your/formats/psrfits.py)
- `your`'s PSRFITS writer sets `fcenter = fch1 + nchans*foff/2` — [your formats/fitswriter.py](https://github.com/thepetabyteproject/your/blob/main/your/formats/fitswriter.py)

### Inferences
- Combining the spec's OBSFREQ definition (centre of channel OBSNCHAN/2, 1-indexed) with `fch1` as the centre of channel 1 gives OBSFREQ = fch1 + (N/2 − 1)·foff. The true band centre is fch1 + (N−1)/2·foff, so the two differ by half a channel. `your`'s `fch1 + N·foff/2` matches neither exactly (it is the band centre + foff/2), so any PSRFITS↔.fil round trip can introduce a ±½-channel frequency shift. DAT_FREQ should be treated as authoritative over OBSFREQ/OBSBW.
- Because the spec says search data "are normally analysed directly without application of the scale factors and offsets", whether a tool applies DAT_SCL/DAT_OFFS changes the bandpass shape of the output. Converting PSRFITS→filterbank with scales applied and then requantising to 8 bits can clip. Converting without them discards the bandpass information permanently, since filterbank has no slot for it.

### Gaps
- I did not read the Hotan et al. (2004) PASA paper itself (ADS: https://ui.adsabs.harvard.edu/abs/2004PASA...21..302H). All keyword detail comes from the v6.1 online documentation.
- I did not inspect DSPSR source directly. Its ZERO_OFF/2-bit behaviour is cited only through sigpyproc3's comments that mirror it.

---

## 3. PSRDADA, GUPPI RAW, VDIF

### Takeaway
DADA files and ring buffers carry a fixed-size ASCII header of newline-separated `KEY value` pairs (size given by `HDR_SIZE`, commonly 4096 B), followed by raw data. Meaning depends on keys like NCHAN/NBIT/NPOL/NDIM/TSAMP/BW/FREQ/UTC_START/OBS_OFFSET/RESOLUTION. GUPPI RAW uses 80-char FITS-like header cards ending in `END`, followed by a data block. VDIF is a framed VLBI format with 32-byte (or legacy 16-byte) headers.

### Cited Findings
- The header is ASCII: "Data attributes are stored as keyword-value pairs, each separated by a new line", accessed with `ascii_header_set(header, keyword, format, ...)` / `ascii_header_get(...)` (fprintf/fscanf-style). "`ascii_header_set` has no idea about the size of the ASCII header", so overflow is the caller's responsibility — [psrdada doc/header.tex](https://sourceforge.net/p/psrdada/code/ci/master/tree/doc/header.tex)
- Standard keys listed in the spec: HEADER, HDR_VERSION, HDR_SIZE, INSTRUMENT, PRIMARY, HOSTNAME, FILE_NAME, FILE_SIZE, FILE_NUMBER, OBS_ID, UTC_START, MJD_START, OBS_OFFSET ("number of bytes from the start of the observation" of the first data byte in this file/stream), OBS_OVERLAP, TELESCOPE, SOURCE, FREQ, BW, NPOL, NBIT, NDIM, TSAMP, RA, DEC — [psrdada doc/header.tex](https://sourceforge.net/p/psrdada/code/ci/master/tree/doc/header.tex)
- Example header shipped with psrdada: `HDR_VERSION 1.0`, `HDR_SIZE 4096`, `BW -400` (MHz; negative = inverted band), `FREQ 1382` (centre MHz), `TELESCOPE PKS`, `INSTRUMENT CASPSR`, `NBIT 8`, `NCHAN 1`, `NDIM 1`, `NPOL 2`, `OBS_OFFSET 0`, `UTC_START 2010-04-13-02:05:45`, `TSAMP 0.00125`, `RESOLUTION 4` — [psrdada examples/src/header.dada](https://sourceforge.net/p/psrdada/code/ci/master/tree/examples/src/header.dada)
- Ring-buffer defaults: `DADA_DEFAULT_BLOCK_KEY 0x0000dada`, `DADA_DEFAULT_BLOCK_NUM 4`, `DADA_DEFAULT_BLOCK_SIZE = page_size*128`, `DADA_DEFAULT_HEADER_SIZE = sysconf(_SC_PAGE_SIZE)` (4096 on typical x86 Linux). The Header Block is a separate ring buffer with the same semaphore interface as the Data Block, and "each sub-block written to the ring buffer represents a unique header" — [psrdada src/dada_def.h](https://sourceforge.net/p/psrdada/code/ci/master/tree/src/dada_def.h), [doc/header.tex](https://sourceforge.net/p/psrdada/code/ci/master/tree/doc/header.tex)
- heimdall's DADA reader needs NCHAN, BANDWIDTH (falls back to BW), CFREQ (falls back to FREQ), NPOL, NBIT, TSAMP, BEAM and UTC_START. It computes `f0 = cfreq − bw/2 + (bw/nchan)/2`, `df = bw/nchan`, so a negative BW yields a descending band — [heimdall Formats/PSRDadaRingBuffer.C](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Formats/PSRDadaRingBuffer.C)
- `your` has a DADA reader/writer module (`your/formats/dada.py`) — [your repo](https://github.com/thepetabyteproject/your/tree/main/your/formats)
- GUPPI RAW: each frame is an ASCII header of 80-character entries that "always ends with 'END' followed by 77 spaces", then a binary block. Files typically hold 128 frames. Samples are 8-bit (Baseband supports only 8-bit). Each channel's data ends with overlap samples repeating the start of the next block. PKTSIZE must give an integer number of packets per payload excluding overlap — [Baseband GUPPI docs](https://baseband.readthedocs.io/en/stable/guppi/index.html)
- VDIF: introduced 2009, spec release 1.1.1 at vlbi.org. Headers are eight 32-bit words (32 B), or four words (16 B) for legacy. The first four words are common, the last four are set by the Extended Data Version (EDV 0–4, plus Mark5B-over-VDIF 0xAB). Headers and payloads are little-endian. Payload is 32 B to about 134 MB. Frames sharing channels form a "thread", and a "frameset" groups all threads covering one time segment — [Baseband VDIF docs](https://baseband.readthedocs.io/en/stable/vdif/index.html)

### Inferences
- The meaning of `RESOLUTION` (minimum byte granularity for splitting the stream) and `ORDER` (e.g. TF/FT/TFP sample ordering, used by MeerKAT/UTMOST-era instruments) is defined by downstream instruments and DSPSR rather than by the core header.tex key list.

### Gaps
- RESOLUTION and ORDER are not defined in psrdada's `doc/header.tex`. I found no primary definition in the sources read.
- The exact GUPPI header keys (BLOCSIZE, OBSNCHAN, NPOL, NBITS, TBIN, PKTIDX, OVERLAP, DIRECTIO) and the VDIF word-level bit fields (seconds from epoch, reference epoch, frame number, frame length in 8-byte units, log2(nchan), bits/sample−1, complex flag, thread/station ID) were not covered by the fetched pages. The VDIF spec PDF at https://vlbi.org/vlbi-standards/vdif/ should be consulted directly.

---

## 4. PRESTO formats (.inf, .dat, .fft, .mask, .rfi/.stats/.bytemask, ACCEL, .singlepulse, .pfd, polycos)

### Takeaway
PRESTO pairs raw native-endian binary arrays (`.dat` float32, `.fft` complex64) with a human-readable `.inf` text sidecar. rfifind writes a compact binary mask (doubles + ints + ragged int lists). accelsearch writes a binary array of `fourierprops` structs (`.cand`) plus a text table. prepfold `.pfd` is a sequential binary dump with ints, length-prefixed strings, doubles and a float64 profile cube.

### Cited Findings

**.inf (text; `writeinf()` in `ioinf.c`, struct `infodata` in `makeinf.h`)** — keys are fixed left-hand labels, `=`, then value:
`Data file name without suffix`; `Telescope used`; `Instrument used`; `Object being observed`; `J2000 Right Ascension (hh:mm:ss.ssss)`; `J2000 Declination     (dd:mm:ss.ssss)`; `Data observed by`; `Epoch of observation (MJD)` (written `%d.%s`, integer + fractional string, so full precision is kept); `Barycentered?           (1 yes, 0 no)`; `Number of bins in the time series` (%-11.0f); `Width of each time series bin (sec)` (%.15g); `Any breaks in the data? (1 yes, 0 no)` with `On/Off bin pair #%3d` lines if yes; `Type of observation (EM band)`; then for Radio: `Beam diameter (arcsec)`, `Dispersion measure (cm-3 pc)`, `Central freq of low channel (MHz)`, `Total bandwidth (MHz)`, `Number of channels`, `Channel bandwidth (MHz)` (X-ray/optical variants exist); `Data analyzed by`; `Any additional notes:` — [presto src/ioinf.c](https://github.com/scottransom/presto/blob/master/src/ioinf.c), [include/makeinf.h](https://github.com/scottransom/presto/blob/master/include/makeinf.h)
- sigpyproc reads both `(1 yes, 0 no)` and an older `(1=yes, 0=no)` barycentric spelling it used to write — [sigpyproc3 io/params.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/io/params.py)

**.fft packing**
- PRESTO's real FFT stores N/2 complex values with the Nyquist term in the imaginary part of bin 0: `/* Set data[0].i to Nyquist value */` in `fastffts.c` and "special cases of freq=0 and Nyquist freq" in `twopass_real_fwd.c` — [presto src/fastffts.c](https://github.com/scottransom/presto/blob/master/src/fastffts.c), [twopass_real_fwd.c](https://github.com/scottransom/presto/blob/master/src/twopass_real_fwd.c)

**.mask (binary, `write_mask()` in `mask.c`)**, in order:

| Field | Type | Meaning |
|---|---|---|
| timesigma | double | time-domain sigma cutoff |
| freqsigma | double | frequency-domain sigma cutoff |
| mjd | double | MJD of time zero |
| dtint | double | seconds per interval |
| lofreq | double | MHz of low channel |
| dfreq | double | channel width MHz |
| numchan | int32 | channels |
| numint | int32 | intervals |
| ptsperint | int32 | samples per interval |
| num_zap_chans | int32 | then `int32[num_zap_chans]` zap_chans (if >0) |
| num_zap_ints | int32 | then `int32[num_zap_ints]` zap_ints (if >0) |
| num_chans_per_int | `int32[numint]` | channels zapped per interval |
| chans[ii] | `int32[num_chans_per_int[ii]]` | written **only if 0 < n < numchan** (fully-zapped intervals have no list) |

Source: [presto src/mask.c](https://github.com/scottransom/presto/blob/master/src/mask.c), [include/mask.h](https://github.com/scottransom/presto/blob/master/include/mask.h)
- `.bytemask` is a raw `numint*numchan` uint8 array with bit flags GOODDATA 0x00, PADDING 0x01, OLDMASK 0x02, USERCHAN 0x04, USERINTS 0x08, BAD_POW 0x10, BAD_STD 0x20, BAD_AVG 0x40 — [presto rfifind.c](https://github.com/scottransom/presto/blob/master/src/rfifind.c), [mask.h](https://github.com/scottransom/presto/blob/master/include/mask.h)
- `.stats` holds int32 numchan, numint, ptsperint, lobin, numbetween, then float32 `datapow[numchan*numint]`, `dataavg[...]` and `datastd[...]`. `.rfi` holds int32 numchan, numint, ptsperint, lobin, numbetween, harmsum, numrfi, float32 fracterror, freqsigma, then numrfi records of {float freq_avg, float freq_var, float sigma_avg, int numobs, bitmask times[ceil(numint/8)], bitmask chans[ceil(numchan/8)]} — [presto rfifind.c](https://github.com/scottransom/presto/blob/master/src/rfifind.c), [rfi_utils.c](https://github.com/scottransom/presto/blob/master/src/rfi_utils.c)

**accelsearch outputs**
- Names: `<root>_ACCEL_<zmax>` (text), `<root>_ACCEL_<zmax>.cand` (binary), `<root>_ACCEL_<zmax>.txtcand`. With jerk search: `_ACCEL_<zmax>_JERK_<wmax>` — [presto accel_utils.c](https://github.com/scottransom/presto/blob/master/src/accel_utils.c)
- The `.cand` file is `fwrite(props, sizeof(fourierprops), numcands)`, a raw array of C structs: double r, float rerr, double z, float zerr, double w, float werr, then floats pow, powerr, sig, rawpow, phs, phserr, cen, cenerr, pur, purerr, locpow — [presto accelsearch.c](https://github.com/scottransom/presto/blob/master/src/accelsearch.c), [include/presto.h](https://github.com/scottransom/presto/blob/master/include/presto.h)
- Text ACCEL fundamentals columns: Cand, Sigma, Summed Power, Coherent Power, Num Harm, Period (ms), Frequency (Hz), FFT 'r' (bin), Freq Deriv (Hz/s), FFT 'z' (bins), FFT 'w' (bins), Accel (m/s²), Notes. Harmonics section: Cand, Harm, Sigma, Power/Loc Pow, Raw Power, FFT 'r', Pred 'r', FFT 'z', Pred 'z', FFT 'w', Pred 'w', Phase (rad), Centroid (0-1), Purity (<p>=1), Notes — [presto accel_utils.c](https://github.com/scottransom/presto/blob/master/src/accel_utils.c)

**.singlepulse**
- Header line `# DM      Sigma      Time (s)     Sample    Downfact` — [presto bin/single_pulse_search.py](https://github.com/scottransom/presto/blob/master/bin/single_pulse_search.py)

**.pfd (prepfold binary; reader `python/presto/prepfold.py`)**, sequential layout:
1. int32 ×5: numdms, numperiods, numpdots, nsub, npart. The reader byte-swaps if any \|value\| > 100000.
2. int32 ×7: proflen, numchan, pstep, pdstep, dmstep, ndmfact, npfact.
3. Four length-prefixed strings (int32 len + bytes): filenm, candnm, telescope, pgdev.
4. char[16] rastr, char[16] decstr.
5. double dt, startT; double ×7: endT, tepoch, bepoch, avgvoverc, lofreq, chan_wid, bestdm.
6. Three groups (topo, bary, fold) of {float pow, float pad; double p1, p2, p3}.
7. double ×7 orbit: p, e, x, w, t, pd, wd.
8. double dms[numdms], periods[numperiods], pdots[numpdots].
9. double profs[npart][nsub][proflen].
10. foldstats[npart][nsub], each 7 doubles: numdata, data_avg, data_var, numprof, prof_avg, prof_var, redchi.

The Python reader also applies GBT SPIGOT/GUPPI-era lofreq/epoch corrections when it reads. Sources: [presto python/presto/prepfold.py](https://github.com/scottransom/presto/blob/master/python/presto/prepfold.py), [include/presto.h](https://github.com/scottransom/presto/blob/master/include/presto.h)

**polycos**
- PRESTO's `polycos.py` parses TEMPO1-format polyco entries (default span 60 min) and notes that polyco generation now uses tempo2 observatory names — [presto python/presto/polycos.py](https://github.com/scottransom/presto/blob/master/python/presto/polycos.py)

### Inferences
- `.dat`, `.fft` and `.cand` have no header or endianness marker. They are native-endian (little-endian on all current platforms), and `.cand` depends on C struct padding. On LP64 with natural alignment, `fourierprops` is 88 bytes (8 + 4 + 4 pad + 8 + 4 + 4 pad + 8 + 12×4). A NumPy dtype must include the pad bytes.
- Because `.mask` omits per-interval lists when an interval is fully zapped, a reader must reproduce that conditional exactly or it will desynchronise.

### Gaps
- I did not verify the `.dat` float32 claim in source this session (`prepdata.c` writes via `chkfwrite` of float arrays). It is consistent with sigpyproc's `.inf`/`.dat` support but not quoted here.
- The `.txtcand` format was not inspected.

---

## 5. heimdall .cand, FETCH input (HDF5 / CSV), Arrow/Parquet

### Takeaway
heimdall writes one tab-separated `.cand` text file per gulp with 9 columns. `your_candmaker.py` turns a CSV candidate list into per-candidate HDF5 files containing `data_freq_time` and `data_dm_time` (default 256×256), which is the input format FETCH consumes. None of the inspected codebases use Arrow/Parquet.

### Cited Findings
- heimdall `.cand` filename: `<output_dir>/<UTC timestamp>_<beam>.cand`. Each line holds `snr \t samp_idx \t samp_idx*dt (s) \t filter_index \t dm_index \t dm \t members \t begin_sample \t end_sample`. Here samp_idx = first_idx + group peak index, filter width = 2^filter_index samples, and members = number of clustered events — [heimdall Pipeline/pipeline.cu](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Pipeline/pipeline.cu), [label_candidate_clusters.cu](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Pipeline/label_candidate_clusters.cu) (`filter_width = 1 << ...`)
- sigpyproc3's candidate table knows heimdall's extra columns `("filter_index", "begin_bin", "end_bin", "members")`. Its pulse CSV export has dm, snr, time_s, bin, width_bins and width_s, and notes that the bin convention is "center for PRESTO and Heimdall, start for matched filter" — [sigpyproc3 core/candidates.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/core/candidates.py)
- `your_candmaker.py` input CSV columns: `file`, `snr`, `width` (**log2 of samples**; code uses `2 ** row["width"]`), `dm`, `label`, `stime` (s), `chan_mask_path`, `num_files`. Defaults are `--frequency_size 256` and `--time_size 256`. Time is decimated by `width//2` and cropped to centre the pulse — [your bin/your_candmaker.py](https://github.com/thepetabyteproject/your/blob/main/bin/your_candmaker.py)
- Candidate HDF5 (`Candidate.save_h5`) has attrs `cand_id`, `tcand`, `dm`, `dm_opt`, `snr`, `snr_opt`, `width`, `label`, `rfi_mask`, `kill_mask`, `filelist`, all `your_header` fields, `tsamp`, `nchans`, `foff` and `nspectra`. Datasets are `data_freq_time` (dims labelled time, frequency; gzip level 9) and optional `data_dm_time` (dims dm, time) — [your your/candidate.py](https://github.com/thepetabyteproject/your/blob/main/your/candidate.py)
- `your` also wraps heimdall (`your_heimdall.py`, `your/utils/heimdall.py` with `generate_dm_list`) — [your utils/heimdall.py](https://github.com/thepetabyteproject/your/blob/main/your/utils/heimdall.py)
- A full-text search of PRESTO, sigproc, sigpyproc3, your, heimdall and psrdada for `parquet|pyarrow|feather` matched only PRESTO's `pixi.lock` dependency lock file. None of these tools read or write Arrow/Parquet — (search of the source snapshots listed at the top)

### Inferences
- The `width` field means different things across tools: heimdall uses a filter index (log2 of boxcar width), your_candmaker CSV uses log2 width, sigpyproc uses `width_bins`, and PRESTO `.singlepulse` uses Downfact in samples. Converters must normalise this.

### Gaps
- I did not inspect the FETCH repository itself (https://github.com/devanshkv/fetch), so the claim that it reads exactly `data_freq_time` and `data_dm_time` at 256×256 rests on `your`'s candmaker defaults and its documented FETCH pairing, not on FETCH source.

---

## 6. Header fields each tool uses, and conversion pitfalls

### Takeaway
Each tool reads only a subset of header fields, and they disagree on band direction, sub-byte bit order, ZERO_OFF handling, scale/offset/weight application, DAT_FREQ uniformity and OFFS_SUB gaps. Filterbank cannot represent per-channel scales, offsets, weights, non-uniform frequencies or row gaps, so converting PSRFITS to .fil loses that information.

### Cited Findings (fields used)
- **PRESTO / PSRFITS**: STT_IMJD/SMJD/OFFS, CHAN_DM, TBIN, NPOL, POL_TYPE, NSBLK, NSUBOFFS, ZERO_OFF (absolute value), OFFS_SUB, DAT_FREQ, DAT_WTS/OFFS/SCL (applied only if non-trivial), TELESCOP+BACKEND (for the Parkes FB_1BIT bit flip) — [presto psrfits.c](https://github.com/scottransom/presto/blob/master/src/psrfits.c)
- **PRESTO / SIGPROC**: nbits, signed, nchans, nifs, tsamp, fch1, foff, az_start, za_start, ibeam, telescope_id, machine_id, tstart, src_raj, src_dej. N comes from the file size — [presto sigproc_fb.c](https://github.com/scottransom/presto/blob/master/src/sigproc_fb.c)
- **heimdall**: SIGPROC nchans, nbits, ibeam, tsamp, fch1, foff, tstart. DADA NCHAN, BANDWIDTH/BW, CFREQ/FREQ, NPOL, NBIT, TSAMP, BEAM, UTC_START — [heimdall SigprocFile.C](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Formats/SigprocFile.C), [PSRDadaRingBuffer.C](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Formats/PSRDadaRingBuffer.C)
- **sigpyproc3 / PSRFITS**: STT_*, RA/DEC, DAT_FREQ (must be uniform), TBIN, NBITS, NPOL/POL_TYPE, ZERO_OFF (dspsr rules), DAT_SCL/OFFS/WTS, OFFS_SUB (gaps rejected), TSUBINT, TEL_AZ/TEL_ZEN, CHAN_DM, IBEAM, NSUBOFFS (kept, not added to tstart) — [sigpyproc3 io/pfits.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/io/pfits.py)
- **your / PSRFITS**: DAT_SCL, DAT_OFFS, DAT_WTS (enabled if non-trivial, warns if columns change between files), POL_TYPE. **No ZERO_OFF** — [your psrfits.py](https://github.com/thepetabyteproject/your/blob/main/your/formats/psrfits.py)

### Cited Findings (pitfalls)
1. **Sub-byte bit order**: SIGPROC is LSB-first ([pack_unpack.c](https://github.com/SixByNine/sigproc/blob/master/src/pack_unpack.c)) and PSRFITS is MSB-first ("earlier samples in higher-order bits", [ATNF doc](https://www.atnf.csiro.au/research/pulsar/psrfits_definition/PsrfitsDocumentation.html)). Copying packed bytes from PSRFITS into a .fil without repacking scrambles channels within each byte. sigpyproc encodes this with different `bitorder` defaults ([bits.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/io/bits.py)).
2. **Parkes FB_1BIT exception**: PRESTO reverses 8-byte groups for this backend only ([psrfits.c](https://github.com/scottransom/presto/blob/master/src/psrfits.c)).
3. **ZERO_OFF**: spec says subtract it, nominally 2^(NBIT−1)−0.5 ([ATNF](https://www.atnf.csiro.au/research/pulsar/psrfits_definition/PsrfitsDocumentation.html)). PRESTO subtracts |ZERO_OFF| for all nbits ([psrfits.c](https://github.com/scottransom/presto/blob/master/src/psrfits.c)). sigpyproc/dspsr skip it for 2- and 16-bit ([pfits.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/io/pfits.py)). `your` never subtracts it ([psrfits.py](https://github.com/thepetabyteproject/your/blob/main/your/formats/psrfits.py)).
4. **Scales/offsets/weights lost**: filterbank has no per-channel scale, offset or weight fields ([read_header.c keyword list](https://github.com/SixByNine/sigproc/blob/master/src/read_header.c)). Per-row DAT_SCL/DAT_OFFS changes (they vary per SUBINT row per the column definitions in the [ATNF doc](https://www.atnf.csiro.au/research/pulsar/psrfits_definition/PsrfitsDocumentation.html)) must be baked into the samples or dropped.
5. **Non-uniform DAT_FREQ**: SIGPROC supports `fchannel` tables only up to 4096 channels, and most readers (heimdall, PRESTO) use only fch1/foff ([header.h](https://github.com/SixByNine/sigproc/blob/master/src/header.h), [SigprocFile.C](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Formats/SigprocFile.C)). PRESTO only warns on non-uniform spacing, and sigpyproc rejects it ([psrfits.c](https://github.com/scottransom/presto/blob/master/src/psrfits.c), [pfits.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/io/pfits.py)).
6. **Band direction**: PSRFITS CHAN_BW<0 means descending frequency ([ATNF](https://www.atnf.csiro.au/research/pulsar/psrfits_definition/PsrfitsDocumentation.html)). PRESTO flips PSRFITS when hi<lo and always reverses .fil channels assuming foff<0 ([sigproc_fb.c](https://github.com/scottransom/presto/blob/master/src/sigproc_fb.c)). sigpyproc keeps DAT_FREQ order with no reversal ([pfits.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/io/pfits.py)). The SIGPROC writer handles inversion by rewriting fch1/foff ([filterbank_header.c](https://github.com/SixByNine/sigproc/blob/master/src/filterbank_header.c)).
7. **Centre-frequency convention**: OBSFREQ is defined as the centre of channel OBSNCHAN/2, not the band centre ([ATNF](https://www.atnf.csiro.au/research/pulsar/psrfits_definition/PsrfitsDocumentation.html)). `your` writes `fch1 + nchans*foff/2` ([fitswriter.py](https://github.com/thepetabyteproject/your/blob/main/your/formats/fitswriter.py)). heimdall derives f0 from CFREQ − BW/2 + chan/2 ([PSRDadaRingBuffer.C](https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Formats/PSRDadaRingBuffer.C)).
8. **Time precision**: PSRFITS splits the start into int MJD + int s + double fractional s ([ATNF](https://www.atnf.csiro.au/research/pulsar/psrfits_definition/PsrfitsDocumentation.html)), while SIGPROC `tstart` is a single double MJD ([read_header.c](https://github.com/SixByNine/sigproc/blob/master/src/read_header.c)). PRESTO `.inf` writes MJD as int + fractional string to keep precision ([ioinf.c](https://github.com/scottransom/presto/blob/master/src/ioinf.c)).
9. **OFFS_SUB gaps / NSUBOFFS**: OFFS_SUB is the row *centre* time ([ATNF](https://www.atnf.csiro.au/research/pulsar/psrfits_definition/PsrfitsDocumentation.html)). PRESTO prefers OFFS_SUB over NSUBOFFS to detect dropped rows and pads them ([psrfits.c](https://github.com/scottransom/presto/blob/master/src/psrfits.c)). sigpyproc refuses gapped files ([pfits.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/io/pfits.py)). A .fil has no way to mark gaps.
10. **NSTOT / partial last row**: the last SUBINT row may be partly filled, and NSTOT gives valid samples per file ([ATNF](https://www.atnf.csiro.au/research/pulsar/psrfits_definition/PsrfitsDocumentation.html)). Converters that write whole rows add junk samples.
11. **Polarisation summing**: PRESTO sums AA+BB without scaling ([psrfits.c](https://github.com/scottransom/presto/blob/master/src/psrfits.c)). `your` averages them (×0.5) ([psrfits.py](https://github.com/thepetabyteproject/your/blob/main/your/formats/psrfits.py)). sigpyproc forms P0+P1 ([pfits.py](https://github.com/FRBs/sigpyproc3/blob/main/src/sigpyproc/io/pfits.py)). Absolute levels therefore differ by 2× between tools.
12. **send_coords always writes coordinates**, including sentinel 0/−1 values ([send_stuff.c](https://github.com/SixByNine/sigproc/blob/master/src/send_stuff.c)). The `nsamples` keyword is ignored on read ([read_header.c](https://github.com/SixByNine/sigproc/blob/master/src/read_header.c)), so the file size is authoritative.

### Inferences
- A safe PSRFITS→.fil converter should: unpack MSB-first; subtract ZERO_OFF following a declared policy; decide explicitly whether to apply DAT_SCL/DAT_OFFS/DAT_WTS (and record that choice somewhere, e.g. in `rawdatafile` or a sidecar); verify uniform DAT_FREQ; fill or flag OFFS_SUB gaps; truncate to NSTOT; write `fch1` = DAT_FREQ[0] and `foff` = DAT_FREQ[1]−DAT_FREQ[0] (preferably reordered so that foff<0 for PRESTO); repack LSB-first for nbits<8; and compute `tstart = STT_IMJD + (STT_SMJD + STT_OFFS)/86400` (+ NSUBOFFS·NSBLK·TBIN/86400 for later split files, after confirming the writer's convention).

### Gaps
- DSPSR's own readers (`dspsr/Kernel/Formats/fits`, `sigproc`) were not inspected directly. Its behaviour is inferred through sigpyproc3's "dspsr unpacker rules" comments.
- Whether STT_* in later time-split PSRFITS files already includes the NSUBOFFS offset varies by writer. sigpyproc explicitly does not add it, and I found no authoritative ruling.
