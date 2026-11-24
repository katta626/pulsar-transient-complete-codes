import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
sys.path.insert(0, str(ROOT / "python"))
from test_sp_cli import write_fil
import filterbank as fb

FB = str(ROOT / "target" / "release" / "fb")
d = Path(tempfile.mkdtemp())
rng = np.random.default_rng(2)
nchan, nsamp, tsamp, P, dm = 32, 200000, 1e-3, 0.4567, 20.0
fr = 1500 - 4.0 * np.arange(nchan)
t = np.arange(nsamp) * tsamp
x = rng.normal(100, 8, (nsamp, nchan))
for c in range(nchan):
    ph = ((t - fb.K_DM * dm * (fr[c] ** -2 - fr[0] ** -2)) / P) % 1.0
    x[:, c] += 3 * np.exp(-0.5 * ((ph - 0.3) / 0.01) ** 2)
write_fil(d / "p.fil", x, 1500.0, -4.0, tsamp)
subprocess.run([FB, "fold", d / "p.fil", "--period", str(P), "--dm", str(dm), "--nsub", "8", "--nbins", "128", "-o", d / "p.fold.json"], check=True, capture_output=True)
subprocess.run([FB, "toa", d / "p.fold.json", "--site", "gmrt", "-o", d / "p.tim"], check=True, capture_output=True)
lines = open(d / "p.tim").read().split("\n")[1:-1]
mjds = np.array([float(l.split()[2]) for l in lines])
errs = np.array([float(l.split()[3]) for l in lines])
turns = (mjds - 60000.0) * 86400 / P
resid = (turns - np.round(turns)) * P * 1e6
print("TOA residuals (us):", np.round(resid, 1), " median error (us):", np.median(errs))
ok = len(lines) == 8 and np.all(np.abs(resid - resid.mean()) < 4 * errs) and np.std(resid) < 3 * np.median(errs)
print("ok" if ok else "FAIL")
sys.exit(0 if ok else 1)
