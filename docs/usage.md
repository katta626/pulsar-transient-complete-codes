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

If `nvcc` is found, it also builds `cuda/build/libfbcuda.so`.

To use the Python package from anywhere, either `pip install -e filterbank` or set `PYTHONPATH=filterbank/python`.

## Command line

| Command | What it does |
|---|---|
| `fb header FILE` | Print the canonical header (JSON): exact start time, the frequency of every channel, scaling policy, any data gaps |
| `fb convert IN... -o OUT.fil` | PSRFITS/DADA/filterbank → SIGPROC filterbank, with correct bit order, scales/offsets/weights, gap filling and summed polarisations |
| `fb rfi IN... -o mask.json [--presto out.mask]` | Build an RFI mask (IQRM, `--sk`, `--stats`, `--zero-dm`); can also write a PRESTO `.mask` |
| `fb export-mask DATA MASK -o OUT` | Convert masks between JSON and PRESTO `.mask` for a given data file |
| `fb sp IN... --dm 0:1000` | Single-pulse search: RFI → dedispersion (`--dedisp auto/brute/subband/fdmt`) → boxcar → clustering. Writes `.cands.csv` and `.singlepulse` (and `.cand` with `--format heimdall`) |
| `fb periodic IN... --dm 0:100 [--accel 50] [--ffa]` | Periodicity search: RFI → batched dedispersion → FFT/acceleration and/or FFA → sifting → folding of the top candidates. Writes `.periodic.csv` and `.folds.json` |
| `fb fft SERIES` / `fb ffa SERIES` | Search one dedispersed series: a PRESTO `.dat/.inf` pair or a SIGPROC `.tim` |
| `fb fold FILE --period P [--pdot] [--dm]` | Fold a series or a filterbank; writes the profile, sub-integrations and reduced χ² |
| `fb toa FOLD.json [--template T.json] --site gmrt` | TOAs per sub-integration by Fourier template matching → tempo2 `.tim` |

Every command accepts PSRFITS (one file or several, sorted by start time), SIGPROC `.fil` and PSRDADA `.dada`. The format is detected from each file's first bytes.

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
- **Backends:** `fb.backend.available()` returns `['cpu']` or `['cpu', 'cuda']`.

## CUDA

On a Linux machine with an NVIDIA GPU:

```bash
cd filterbank && make cuda
```

- This builds `cuda/build/libfbcuda.so`, or set `FB_CUDA_LIB` to its path.
- `fb.dedisp.dedisperse(..., backend="cuda")` then runs dedispersion on the GPU, using the same delay table as the CPU code, so both give identical integer delays.
- Zero-DM and boxcar kernels are also available, as `fb.cuda.zero_dm` and `fb.cuda.boxcar_snr`.

## Where things are

- [reference/](reference/00-overview.md): what each upstream package does, the formats, and the math.
- [reference/06-our-plan.md](reference/06-our-plan.md): the module plan.
- [CREDITS.md](CREDITS.md): where each borrowed idea or format came from.
- [math-maps/](math-maps/README.md): one map per library of the equations it uses and where they are in its code.
