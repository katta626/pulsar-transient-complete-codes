import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests" / "cuda_emu"))
import build as emu

os.environ["FB_CUDA_LIB"] = str(emu.build(tempfile.mkdtemp()))
sys.path.insert(0, str(ROOT / "python"))
import filterbank as fb
from filterbank import cuda as CU

fails = 0


def check(name, cond):
    global fails
    print(("ok   " if cond else "FAIL ") + name)
    fails += 0 if cond else 1


check("emulated device visible", fb.backend.available() == ["cpu", "cuda"])

rng = np.random.default_rng(1)
nchan, nsamp, tsamp = 32, 3000, 1e-3
freqs = 1500.0 - 2.0 * np.arange(nchan)
h = fb.Header({
    "source_name": "EMU", "telescope": "Fake", "backend": "FAKE", "ra_deg": 0.0, "dec_deg": 0.0, "az_deg": None, "za_deg": None,
    "start": {"imjd": 60000, "sec": 0.0}, "tsamp": tsamp, "nsamples": nsamp, "freqs": {"mhz": list(freqs), "chan_bw": -2.0},
    "npol": 1, "pol_type": "Intensity", "nbits": 32, "scaling": {"zero_off": 0.0, "signed": False, "applied": True},
    "barycentric": False, "beam": 0, "nbeams": 1, "order": "TimePolFreq", "provenance": [],
})
data = rng.normal(0, 1, (nsamp, nchan)).astype(np.float32)
kill = np.zeros(nchan, dtype=bool)
kill[5] = True
p = fb.dedisp.plan(h, 0, 150)
cpu = fb.dedisp.dedisperse(h, p, data, killmask=kill, backend="cpu", method="brute")
gpu = fb.dedisp.dedisperse(h, p, data, killmask=kill, backend="cuda")
check("dedisperse kernel matches CPU brute force", cpu.shape == gpu.shape and np.allclose(cpu, gpu, atol=1e-3))

z = CU.zero_dm(data)
check("zero-DM kernel matches numpy", np.allclose(z, fb.rfi.zero_dm(data), atol=1e-4))

series = rng.normal(0, 1, (4, 3000)).astype(np.float32)
widths = np.array([1, 2, 4, 8, 16], dtype=np.int32)
snr, bw = CU.boxcar_snr(series, widths)
cs = np.concatenate([np.zeros((4, 1)), np.cumsum(series, axis=1, dtype=np.float64)], axis=1)
ref = np.full(series.shape, -1e30)
for w in widths:
    s = cs[:, w:] - cs[:, :-w]
    sub = s[:, :: max(1, s.shape[1] // 16384)]
    mid = sub.shape[1] // 2
    med = np.partition(sub, mid, axis=1)[:, mid : mid + 1]
    sd = 1.4826 * np.partition(np.abs(sub - med), mid, axis=1)[:, mid : mid + 1]
    s = (s - med) / sd
    ref[:, : s.shape[1]] = np.maximum(ref[:, : s.shape[1]], s)
check("boxcar kernel matches per-width noise reference", np.allclose(snr, ref, atol=2e-3))

fd_cpu = fb.dedisp.dedisperse(h, p, data, killmask=kill, backend="cpu", method="fdmt")
fd_gpu = fb.dedisp.dedisperse(h, p, data, killmask=kill, backend="cuda", method="fdmt")
check("FDMT kernels match CPU FDMT", fd_cpu.shape == fd_gpu.shape and np.allclose(fd_cpu, fd_gpu, atol=1e-3))

n, ts, f0 = 1 << 14, 2e-3, 23.7
tt = np.arange(n) * ts
ph = (f0 * tt) % 1.0
rows = rng.normal(0, 1, (3, n)).astype(np.float32)
rows[1] += 0.6 * (np.abs(ph - 0.5) < 0.04)
res = CU.fft_search(rows, ts, f_min=1.0, f_max=200.0, max_harm=8, sigma=6.0)
top = res[1][0] if res[1] else None
check("GPU FFT search finds the pulsar row", top is not None and abs(top["f_hz"] - f0) < 1.0 / (n * ts) and top["sigma"] > 8)
check("GPU FFT search quiet on noise rows", all(not r or r[0]["sigma"] < 7 for r in (res[0], res[2])))

sys.exit(1 if fails else 0)
