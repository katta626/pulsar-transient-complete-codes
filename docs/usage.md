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
  - The spectrum is whitened with a window of up to `8·zmax` bins, so a strongly drifting signal does not normalise itself away.
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
  - It uses the psrdada client code, so it should work alongside the real tools (`dada_db`, `dada_dbdisk`, ...). That has not been tested yet.
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
- The library compiles with `nvcc` 12.0 (in the `cuda` Linux system, see `~/docker/README.md`). It has not been run on a GPU.
- **Without a GPU:** `tests/test_cuda_emu.py` compiles `fbcuda.cu` against small stand-in CUDA headers and runs the kernels on CPU threads. It compares each one with the CPU code, so kernel logic can be checked on a Mac. It is slow, about 30 s for the tests, and is for correctness only.

## Benchmarks

```bash
python filterbank/bench/bench.py --nchan 512 --nsamp 524288 --dm-hi 600
```

This writes a synthetic 8-bit filterbank with a dispersed burst and a pulsar, or uses `--file obs.fil`. It then times each step and reports what was found:

- `fb rfi`, `fb sp` (subband and FDMT) and `fb periodic`;
- PRESTO (`rfifind`, `prepsubband`, `single_pulse_search.py`, `realfft` + `accelsearch`) and heimdall, but only if they are on `PATH`.

The table goes to the terminal and the results to `bench_work/results.json`. The PRESTO steps use `single_pulse_search.py -b` and `rednoise` (see the comparison below).

On an M1 Pro, CPU only (512 channels, 33.5 s of data, DM 0–600):

| step | seconds | found |
|---|---|---|
| `fb rfi` | 0.7 | |
| `fb sp --dedisp subband` | 13.9 | DM 349.8, S/N 52 at 10.000 s (injected DM 350, 10 s) |
| `fb sp --dedisp fdmt` | 21.9 | DM 349.8, S/N 53 |
| `fb periodic --dm 0:120` | 9.2 | P = 0.2347 s at DM 58 (injected 0.2347 s, DM 60) |

## Compared with PRESTO

Synthetic data, same input files for both. PRESTO 6.3.3 (commit `7f66332`, gcc -mcpu=native, OpenMP) in an arm64 Linux container on an M1 Pro (6 fast + 2 slow cores); `fb` runs directly on macOS, so the container's overhead counts against PRESTO. Times are wall-clock, best of 3 (two runs for `single_pulse_search.py`, one for the `nsub 32` run). FFTW planning is excluded for PRESTO: without saved plans it takes about 54 s here, because FFTW lacks a fine timer in this ARM build.

| Step | PRESTO 1 thread | PRESTO 8 threads | fb 1 thread | fb 8 threads |
|---|---|---|---|---|
| FFT search, no acceleration, 8 harmonics | 0.01 s | 0.01 s | 0.03 s | 0.02 s |
| Acceleration search, zmax 200, 8 harmonics | 0.79 s | 0.15 s | 1.25 s | 0.28 s |
| Jerk search, zmax 100, wmax 300, 4 harmonics | 7.1 s | 1.2 s | 4.5 s | 1.0 s |
| RFI mask, 512 channels, 33.5 s | 6.2 s | 1.9 s | 0.21 s | 0.19 s |

The `fb` acceleration and jerk commands include reading the series, whitening and the forward FFT; `accelsearch` starts from a finished `.fft`. The RFI methods differ (`rfifind` statistics against IQRM), so that row does not compare equal work.

Single pulse, DM 0-600, 512 channels, 33.5 s, 8 threads. The DM grids differ (PRESTO 601 trials at step 1, ours 353), and `single_pulse_search.py` is single-threaded:

| | seconds | burst S/N (DM 350, 10 s) |
|---|---|---|
| `prepsubband` (nsub 128) + `single_pulse_search.py -b` | 28 + 122 | 32.7 |
| `prepsubband`, nsub 32, 1 thread (dedispersion only) | 9 | |
| `fb sp` subband, RFI included | 14.0 | 52.0 |
| `fb sp` fdmt, RFI included | 24.9 | 52.7 |

A matched boxcar on PRESTO's own series gives 53.6. Its lower value comes from the burst inflating the block's estimated standard deviation. An older PRESTO (2017) run without `-b` did not report the burst, because bad-block rejection discarded its block. With that version, no method (brute, subband, fdmt, PRESTO) reported anything above 6 on one 33.5 s noise file, DM 0-300.

Sensitivity (sigma after trials; 2^19 white-noise samples plus a pulsar drifting by z bins; 8 harmonics, `zmax 200`, matched frequency ranges; PRESTO's `-flo/-fhi` refer to the highest harmonic):

| z | `accelsearch` | `fb fft --zmax` |
|---|---|---|
| 0 | 21.8 | 22.4 |
| 40 | 18.9 | 20.1 |
| 120 | 9.0 | 9.1 |
| 180 | 9.7 | 9.4 |
| 120, with red noise | 8.0 | 9.8 |

These are single realisations. On five white-noise files PRESTO listed a chance candidate of 2.2 to 3.2 sigma in four; ours listed none above 2 at `zmax 200` (our trials count is stricter).

Not compared: RFI mask quality, folding, jerk sensitivity, heimdall, real telescope data.

## Where things are

- [reference/](reference/00-overview.md): what each upstream package does, the formats, and the math.
- [reference/06-our-plan.md](reference/06-our-plan.md): the module plan.
- [CREDITS.md](CREDITS.md): where each borrowed idea or format came from.
- [math-maps/](math-maps/README.md): one map per library of the equations it uses and where they are in its code.
