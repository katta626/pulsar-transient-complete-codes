import csv
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "tests"))
import filterbank as fb
from test_sp_cli import write_fil

FB = ROOT / "target" / "release" / "fb"
fails = 0


def check(name, cond):
    global fails
    print(("ok   " if cond else "FAIL ") + name)
    fails += 0 if cond else 1


def run(*a):
    r = subprocess.run([str(FB), *map(str, a)], capture_output=True, text=True)
    if r.returncode:
        print(r.stderr)
    return r


d = Path(tempfile.mkdtemp())
rng = np.random.default_rng(4)
nchan, nsamp, tsamp, fch1, foff = 64, 300_000, 1e-3, 1500.0, -4.0
freqs = fch1 + foff * np.arange(nchan)
period, dm = 0.7123, 35.0
t = np.arange(nsamp) * tsamp
data = rng.normal(100, 8, size=(nsamp, nchan))
for c in range(nchan):
    tc = t - fb.K_DM * dm * (freqs[c] ** -2 - freqs[0] ** -2)
    ph = (tc / period) % 1.0
    data[:, c] += 2.5 * (np.abs(ph - 0.5) < 0.015)
write_fil(d / "psr.fil", data, fch1, foff, tsamp)

t0 = time.time()
r = run("periodic", d / "psr.fil", "--dm", "0:70", "--dm-step", "5", "--ffa", "--pmin", "0.3", "--pmax", "3", "--fold", "5", "-o", d / "psr")
print(r.stdout.strip(), f"\n({time.time() - t0:.1f} s wall)")
rows = list(csv.DictReader(open(d / "psr.periodic.csv")))
best = max(rows, key=lambda x: float(x["sigma"]))
p = float(best["period_s"])
harm = min(abs(p / period - k) for k in (0.5, 1, 2, 1 / 3, 3)) < 2e-3
check("periodic: top candidate is the pulsar (or harmonic)", harm and abs(float(best["dm"]) - dm) <= 10)

with fb.open(d / "psr.fil") as rd:
    h = rd.header
    x = rd.read()[:, 0, :]
pl = fb.dedisp.plan(h, dms=[dm])
series = fb.dedisp.dedisperse(h, pl, x)[0]
fb.write_dat(d / "psr_DM35.00", h, dm, series)
r = run("fft", d / "psr_DM35.00.dat", "--fmin", "0.5", "--fmax", "50", "-o", d / "f.csv")
rows = list(csv.DictReader(open(d / "f.csv")))
check("fft on PRESTO .dat finds the period", any(abs(float(x["period_s"]) / period - 1) < 1e-3 for x in rows[:3]))
r = run("ffa", d / "psr_DM35.00.dat", "--pmin", "0.5", "--pmax", "2", "-o", d / "a.csv")
rows = list(csv.DictReader(open(d / "a.csv")))
check("ffa on PRESTO .dat finds the period", abs(float(rows[0]["period_s"]) / period - 1) < 1e-3)
r = run("fold", d / "psr.fil", "--period", period, "--dm", dm, "-o", d / "fold.json")
chi_good = float(r.stdout.split("=")[1].split()[0])
r = run("fold", d / "psr.fil", "--period", period * 1.01, "--dm", dm, "-o", d / "fold2.json")
chi_bad = float(r.stdout.split("=")[1].split()[0])
check(f"fold filterbank: chi2 {chi_good:.1f} on-period vs {chi_bad:.1f} off", chi_good > 5 * chi_bad)
print(f"{fails} failure(s)")
sys.exit(1 if fails else 0)
