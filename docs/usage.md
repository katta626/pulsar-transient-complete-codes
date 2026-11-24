# Using `filterbank`

Everything runs from one Rust build. On a Mac this is CPU only. On Linux with `nvcc`, the CUDA library is built too, and `backend="auto"` uses it.

## Build

Requirements: a Rust toolchain (`cargo`). Python needs only `numpy`. Building the CUDA backend also needs `nvcc`. The tests use `astropy` to write PSRFITS fixtures.

```bash
cd filterbank && make
```

```bash
cd filterbank && make test
```

`make` produces two things:

- `target/release/fb`, the command-line tool;
- `target/release/libfbcore.{dylib,so}`, the shared library used by Python and by C/C++/CUDA through `include/fbcore.h`.

It also builds `dada/build/libfbdada.{dylib,so}`, the PSRDADA ring-buffer client. It is plain C from the vendored psrdada sources and needs no psrdada install. If `nvcc` is found, it also builds `cuda/build/libfbcuda.so`.

To use the Python package from anywhere, either `pip install -e filterbank` or set `PYTHONPATH=filterbank/python`.

## Command line

| Command | What it does |
|---|---|
| `fb header FILE` | Print the canonical header (JSON): exact start time, the frequency of every channel, scaling policy, any data gaps |
| `fb convert IN... -o OUT.fil` | PSRFITS/DADA/filterbank → SIGPROC filterbank, with correct bit order, scales/offsets/weights, gap filling and summed polarisations |
| `fb rfi IN... -o mask.json [--presto out.mask]` | Build an RFI mask (IQRM, `--sk`, `--stats`, `--zero-dm`); can also write a PRESTO `.mask` |
| `fb export-mask DATA MASK -o OUT` | Convert masks between JSON and PRESTO `.mask` for a given data file |
| `fb sp IN... --dm 0:1000` | Single-pulse search: RFI → dedispersion (`--dedisp auto/brute/subband/fdmt`) → boxcar → clustering. Writes `.cands.csv` and `.singlepulse` (and `.cand` with `--format heimdall`) |
| `fb periodic IN... --dm 0:100 [--accel 50 \| --zmax 200 [--wmax 600]] [--ffa]` | Periodicity search: RFI → batched dedispersion → FFT search and/or FFA → sifting → folding of the top candidates. Writes `.periodic.csv` and `.folds.json` |
| `fb fft SERIES` / `fb ffa SERIES` | Search one dedispersed series: a PRESTO `.dat/.inf` pair or a SIGPROC `.tim` |
| `fb fold FILE --period P [--pdot] [--dm]` | Fold a series or a filterbank; writes the profile, sub-integrations and reduced χ² |
| `fb toa FOLD.json [--template T.json] --site gmrt` | TOAs per sub-integration by Fourier template matching → tempo2 `.tim` |

Every command accepts PSRFITS (one file or several, sorted by start time), SIGPROC `.fil` and PSRDADA `.dada`. The format is detected from each file's first bytes.

## Choosing a dedisperser

`--dedisp` (CLI) or `method=` (Python):

- **`brute`:** exact. Best for small channel counts.
- **`subband`:** two-stage. Fast, with a small extra smearing.
- **`fdmt`:** the tree algorithm. Its cost depends on the largest DM delay, not on the number of DM trials, so it wins for fine DM grids.
- **`auto`:** uses subband at 256 channels or more, brute force otherwise.

Timings on an M1 Pro:

| Data | brute | subband | fdmt |
|---|---|---|---|
| 4096 channels, 1001 DMs, 2.1 s of data | 4.4 s | 0.45 s | 0.47 s |
| 256 channels, 1001 DMs (0–500), 8.4 s of data | 0.24 s | — | 0.07 s |

## Acceleration and jerk searches

There are two ways to search for binary pulsars. Both run under `fb periodic` and `fb fft`.

- **`--accel MAX`** (m/s²): time-domain resampling, one FFT per acceleration trial. Cheap for small accelerations.
- **`--zmax Z [--wmax W]`**: Fourier-domain search, like PRESTO `accelsearch`.
  - `z` is the frequency drift in Fourier bins over the observation, and `w` is the drift of `z` (jerk).
  - Step sizes are `--dz 2` and `--dw 20` by default.
  - Like PRESTO, `zmax`/`wmax` apply to the highest harmonic summed.
  - The top `--refine 50` candidates are then optimised in `(f, ḟ, f̈)` by interpolating the spectrum, and the CSV reports the optimised values (at mid-observation).
  - Speed: 2²⁰ samples, 0.5–500 Hz, `zmax 200`, 16 harmonics takes about 5 s on an M1 Pro.
  - A jerk search costs roughly `2·wmax/dw + 1` times more.

## Mixing with PRESTO and other tools

The physical-unit types are what make the tools interchangeable:

- **RFI masks:** a mask is stored in MHz × MJD, not as channel indices.
  - `fb rfi obs.fil --presto obs_rfifind.mask` writes a mask PRESTO's `prepsubband -mask` can use.
  - A PRESTO `rfifind` mask can be passed straight to our searches with `fb sp obs.fil --mask obs_rfifind.mask`. It is re-mapped to the data's own channel order.
- **Dedispersed series:**
  - Ours can be written as `.dat/.inf` (`filterbank.write_dat`), for PRESTO `accelsearch`/`prepfold`.
  - PRESTO's `.dat` files can be read by `fb fft`, `fb ffa` and `fb fold`.
- **Candidates:** written as PRESTO `.singlepulse`, heimdall `.cand` and CSV; heimdall `.cand` can also be read back in Python.

## Python

```python
import filterbank as fb

with fb.open("obs.fits") as r:              # PSRFITS / .fil / .dada
    h = r.header                            # exact epoch, per-channel MHz, gaps
    block = r.read(0, 65536)                # (nsamp, npol, nchan) float32

mask = fb.rfi.find(h, block, interval=1.0)  # IQRM by default; sk=True, stats=True, zero_dm=True
clean, n = fb.rfi.apply(h, block, mask)

plan = fb.dedisp.plan(h, 0, 1000)                           # tolerance-spaced DMs
trials = fb.dedisp.dedisperse(h, plan, clean, method="subband", backend="auto")
cands = fb.singlepulse.search(h, plan, trials, threshold=7)

res = fb.singlepulse.run("obs.fil", dm_lo=0, dm_hi=1000, dedisp="fdmt")   # whole pipeline
scored = fb.classify.classify("obs.fil", res["candidates"], save_dir="cuts")  # cut-outs + score
```

- **Classifier:** `fb.classify.classify(..., model="fetch_like.onnx")` uses an ONNX model if `onnxruntime` is installed. Without one it uses the built-in physical score: broadband fraction, zero-DM ratio and the shape of S/N against DM.
- **Live data:** `fb.realtime.run_live("growing.fil", gulp=65536, dm_hi=2000)` follows a file while the backend is still writing it. `fb.realtime.watch_directory(...)` yields new files as they appear.
- **PSRDADA ring buffers:** `fb.realtime.run_live_ring("dada", gulp=65536, dm_hi=2000)` attaches to a running ring (hex key, as for `dada_db -k dada`) as a reader. It parses the ring's ASCII header like a `.dada` file and searches each gulp as it arrives, until end of data.
  - `fb.ring.follow(key, gulp, overlap)` yields `(start_sample, block)` for your own processing.
  - `fb.ring.create/destroy` make and remove a ring.
  - `fb.ring.play("obs.dada", key)` replays a file into a ring, for testing a live setup offline.
  - It works alongside the real psrdada tools (`dada_db`, `dada_dbdisk`, ...), since it is the same code.
  - On macOS the default System V limits are small (4 MB, 8 segments per process), so use a small ring there.
- **Backends:** `fb.backend.available()` returns `['cpu']` or `['cpu', 'cuda']`.

## CUDA

On a Linux machine with an NVIDIA GPU:

```bash
cd filterbank && make cuda
```

- This builds `cuda/build/libfbcuda.so`, or set `FB_CUDA_LIB` to its path.
- `fb.dedisp.dedisperse(..., backend="cuda")` then runs dedispersion on the GPU, using the same delay table as the CPU code, so both give identical integer delays.
- With `method="fdmt"` it runs FDMT on the GPU instead. It executes the same schedule that the Rust code builds.
- Zero-DM and boxcar kernels are also available, as `fb.cuda.zero_dm` and `fb.cuda.boxcar_snr`.
- `fb.cuda.fft_search(series, tsamp, f_min, f_max, max_harm, sigma)` runs a batched FFT periodicity search over many dedispersed series at once (cuFFT, dereddening, interbinning, harmonic sums). It returns one candidate list per series.
- **Without a GPU:** `tests/test_cuda_emu.py` compiles `fbcuda.cu` against small stand-in CUDA headers and runs the kernels on CPU threads. It compares each one with the CPU code, so kernel logic can be checked on a Mac. It is slow, about 30 s for the tests, and is for correctness only.

## Benchmarks

```bash
python filterbank/bench/bench.py --nchan 512 --nsamp 524288 --dm-hi 600
```

This writes a synthetic 8-bit filterbank with a dispersed burst and a pulsar, or uses `--file obs.fil`. It then times each step and reports what was found:

- `fb rfi`, `fb sp` (subband and FDMT) and `fb periodic`;
- PRESTO (`rfifind`, `prepsubband`, `single_pulse_search.py`, `realfft` + `accelsearch`) and heimdall, but only if they are on `PATH`.

The table goes to the terminal and the results to `bench_work/results.json`. On an M1 Pro, CPU only (512 channels, 33.5 s of data, DM 0–600):

| step | seconds | found |
|---|---|---|
| `fb rfi` | 0.7 | |
| `fb sp --dedisp subband` | 13.9 | DM 349.8, S/N 52 at 10.000 s (injected DM 350, 10 s) |
| `fb sp --dedisp fdmt` | 21.9 | DM 349.8, S/N 53 |
| `fb periodic --dm 0:120` | 9.2 | P = 0.2347 s at DM 58 (injected 0.2347 s, DM 60) |

## Where things are

- [reference/](reference/00-overview.md): what each upstream package does, the formats, and the math.
- [reference/06-our-plan.md](reference/06-our-plan.md): the module plan.
- [CREDITS.md](CREDITS.md): where each borrowed idea or format came from.
- [math-maps/](math-maps/README.md): one map per library of the equations it uses and where they are in its code.
