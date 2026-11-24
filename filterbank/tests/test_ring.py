import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
import filterbank as fb

fails = 0


def check(name, cond):
    global fails
    print(("ok   " if cond else "FAIL ") + name)
    fails += 0 if cond else 1


if not fb.ring.available():
    print("skip: libfbdada not built")
    sys.exit(0)

d = Path(tempfile.mkdtemp())
rng = np.random.default_rng(3)
nchan, nsamp, tsamp, fcen, bw = 128, 100_000, 256e-6, 1364.75, -192.0
freqs = fcen - bw / 2 + bw / nchan * (np.arange(nchan) + 0.5)
data = rng.normal(64, 5, size=(nsamp, nchan))
t = np.arange(nsamp) * tsamp
dm, t0 = 150.0, 11.0
for c in range(nchan):
    data[:, c] += 3.5 * np.exp(-0.5 * ((t - t0 - fb.K_DM * dm * (freqs[c] ** -2 - freqs[0] ** -2)) / 1.5e-3) ** 2)
raw = np.clip(np.rint(data), 0, 255).astype(np.uint8)
text = "\n".join([
    "HDR_VERSION 1.0", "HDR_SIZE 4096", "SOURCE RINGTEST", "TELESCOPE Fake", "INSTRUMENT FAKE",
    f"FREQ {fcen}", f"BW {bw}", f"NCHAN {nchan}", "NPOL 1", "NBIT 8", "NDIM 1",
    f"TSAMP {tsamp * 1e6}", "UTC_START 2026-11-25-10:00:00", "OBS_OFFSET 0",
]) + "\n"
dada = d / "obs.dada"
dada.write_bytes(text.encode().ljust(4096, b"\0") + raw.tobytes())

key = format(0xF000 + os.getpid() % 0x0F00, "x")

fb.ring.create(key, nbufs=4, bufsz=nchan * 2048, hdr_nbufs=1, hdr_bufsz=4096)
try:
    play = [sys.executable, "-c", f"import sys; sys.path.insert(0, {str(ROOT / 'python')!r}); import filterbank as fb; fb.ring.play({str(dada)!r}, {key!r}, chunk=50000)"]
    w = subprocess.Popen(play)
    blocks = list(fb.ring.follow(key, gulp=30000, overlap=1000))
    w.wait(timeout=60)
    got = np.concatenate([blocks[0][1][:, 0, :]] + [b[:, 0, :][1000:] for _, b in blocks[1:]])
    check("ring streams every sample in order", got.shape == raw.shape and np.array_equal(got, raw.astype(np.float32)))
    check("blocks overlap by the requested amount", all(b[0] == i * 30000 for i, b in enumerate(blocks)))

    w = subprocess.Popen(play)
    cands = fb.realtime.run_live_ring(key, gulp=40000, dm_lo=0, dm_hi=300, threshold=8, max_width=64)
    w.wait(timeout=60)
    with fb.open(dada) as r:
        h = r.header
    best = max(cands, key=lambda c: c["snr"]) if cands else None
    check("live ring search finds the burst", best is not None and abs(best["dm"] - dm) < 10 and abs(fb.singlepulse.candidate_time(best, h) - t0) < 0.05)
    hdr, _ = fb.ring.layout(text)
    check("ring header matches the file reader", hdr["freqs"] == h["freqs"] and hdr["start"] == h["start"] and abs(hdr["tsamp"] - h["tsamp"]) < 1e-15)
finally:
    fb.ring.destroy(key)

sys.exit(1 if fails else 0)
