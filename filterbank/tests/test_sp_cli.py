import csv
import struct
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FB = ROOT / "target" / "release" / "fb"
KDM = 1 / 2.41e-4


def write_fil(path, data, fch1, foff, tsamp, tstart=60000.0):
    def s(x):
        return struct.pack("<i", len(x)) + x.encode()

    h = s("HEADER_START") + s("source_name") + s("SYNTH")
    for k, v in (("telescope_id", 0), ("machine_id", 0), ("data_type", 1), ("nchans", data.shape[1]), ("nbits", 8), ("nifs", 1), ("ibeam", 1), ("nbeams", 1)):
        h += s(k) + struct.pack("<i", v)
    for k, v in (("fch1", fch1), ("foff", foff), ("tstart", tstart), ("tsamp", tsamp)):
        h += s(k) + struct.pack("<d", v)
    h += s("HEADER_END")
    Path(path).write_bytes(h + np.clip(np.round(data), 0, 255).astype(np.uint8).tobytes())


def main():
    rng = np.random.default_rng(3)
    nchan, nsamp, tsamp, fch1, foff = 256, 200_000, 128e-6, 1460.0, -0.8
    freqs = fch1 + foff * np.arange(nchan)
    data = rng.normal(64, 6, size=(nsamp, nchan))
    dm, t0, w = 220.0, 11.0, 1.5e-3
    delays = KDM * dm * (freqs ** -2 - freqs[0] ** -2)
    t = np.arange(nsamp) * tsamp
    for c in range(nchan):
        data[:, c] += 4.0 * np.exp(-0.5 * ((t - t0 - delays[c]) / w) ** 2)
    data[:, 40] += 40 * np.abs(rng.normal(size=nsamp))
    d = Path(tempfile.mkdtemp())
    write_fil(d / "x.fil", data, fch1, foff, tsamp)
    t1 = time.time()
    r = subprocess.run([str(FB), "sp", str(d / "x.fil"), "--dm", "0:500", "--threshold", "7", "-o", str(d / "x")], capture_output=True, text=True)
    el = time.time() - t1
    print(r.stdout.strip() or r.stderr.strip())
    if r.returncode != 0:
        return 1
    rows = list(csv.DictReader(open(d / "x.cands.csv")))
    best = max(rows, key=lambda x: float(x["snr"]))
    tbest = (float(best["mjd"]) - 60000.0) * 86400
    ok = abs(float(best["dm"]) - dm) < 10 and abs(tbest - t0) < 0.01 and float(best["snr"]) > 15
    print(f"best: dm={float(best['dm']):.1f} t={tbest:.4f}s snr={float(best['snr']):.1f}  ({el:.2f} s wall for {nsamp*tsamp:.1f} s of data)")
    print("ok" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
