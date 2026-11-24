import argparse
import json
import shutil
import struct
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FB = ROOT / "target" / "release" / "fb"
K_DM = 1.0 / 2.41e-4


def synth(path, nchan, nsamp, tsamp, fch1, foff, burst, pulsar, seed=1):
    rng = np.random.default_rng(seed)
    freqs = fch1 + foff * np.arange(nchan)

    def s(x):
        return struct.pack("<i", len(x)) + x.encode()

    h = s("HEADER_START") + s("source_name") + s("BENCH")
    for k, v in (("telescope_id", 0), ("machine_id", 0), ("data_type", 1), ("nchans", nchan), ("nbits", 8), ("nifs", 1), ("ibeam", 1), ("nbeams", 1)):
        h += s(k) + struct.pack("<i", v)
    for k, v in (("fch1", fch1), ("foff", foff), ("tstart", 60000.0), ("tsamp", tsamp)):
        h += s(k) + struct.pack("<d", v)
    h += s("HEADER_END")
    bdm, bt0, bamp = burst
    pdm, period, pamp = pulsar
    bdel = K_DM * bdm * (freqs ** -2 - freqs[0] ** -2)
    pdel = K_DM * pdm * (freqs ** -2 - freqs[0] ** -2)
    chunk = 1 << 16
    with open(path, "wb") as f:
        f.write(h)
        for s0 in range(0, nsamp, chunk):
            n = min(chunk, nsamp - s0)
            t = (s0 + np.arange(n))[:, None] * tsamp
            x = rng.normal(64, 6, (n, nchan)).astype(np.float32)
            x += bamp * np.exp(-0.5 * ((t - bt0 - bdel) / 1e-3) ** 2)
            ph = ((t - pdel) / period) % 1.0
            x += pamp * (np.abs(ph - 0.5) < 0.02)
            f.write(np.clip(np.rint(x), 0, 255).astype(np.uint8).tobytes())


def run(cmd, cwd):
    t = time.perf_counter()
    p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    dt = time.perf_counter() - t
    ok = p.returncode == 0
    if not ok:
        sys.stderr.write(f"[{' '.join(map(str, cmd[:2]))}] failed:\n{p.stderr[-2000:]}\n")
    return dt, ok


def read_singlepulse(paths):
    rows = []
    for p in paths:
        for line in Path(p).read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            v = line.split()
            rows.append((float(v[0]), float(v[1]), float(v[2])))
    return max(rows, key=lambda r: r[1]) if rows else None


def read_heimdall(paths):
    rows = []
    for p in paths:
        for line in Path(p).read_text().splitlines():
            v = line.split()
            if len(v) >= 6:
                rows.append((float(v[5]), float(v[0]), float(v[2])))
    return max(rows, key=lambda r: r[1]) if rows else None


def main():
    ap = argparse.ArgumentParser(description="time the filterbank pipeline against PRESTO and heimdall where installed")
    ap.add_argument("--file", help="existing 8-bit filterbank to use instead of synthetic data")
    ap.add_argument("--work", default="bench_work")
    ap.add_argument("--nchan", type=int, default=512)
    ap.add_argument("--nsamp", type=int, default=1 << 19)
    ap.add_argument("--tsamp", type=float, default=64e-6)
    ap.add_argument("--dm-hi", type=float, default=600.0)
    ap.add_argument("--burst", default="350,10.0,3.0", help="DM,time_s,amplitude")
    ap.add_argument("--pulsar", default="60,0.2347,1.2", help="DM,period_s,amplitude")
    ap.add_argument("--skip", default="", help="comma list of: fb,presto,heimdall")
    a = ap.parse_args()
    work = Path(a.work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    burst = tuple(float(x) for x in a.burst.split(","))
    pulsar = tuple(float(x) for x in a.pulsar.split(","))
    fil = Path(a.file).resolve() if a.file else work / "bench.fil"
    if not a.file and not fil.exists():
        t = time.perf_counter()
        synth(fil, a.nchan, a.nsamp, a.tsamp, 1500.0, -300.0 / a.nchan, burst, pulsar)
        print(f"synthetic data: {fil} ({fil.stat().st_size / 1e6:.0f} MB) in {time.perf_counter() - t:.1f} s")
    skip = set(filter(None, a.skip.split(",")))
    results = []

    def add(tool, step, dt, ok, found=None):
        results.append({"tool": tool, "step": step, "seconds": round(dt, 3), "ok": ok, "found": found})

    if "fb" not in skip and FB.exists():
        d = work / "fb"
        d.mkdir(exist_ok=True)
        dt, ok = run([FB, "rfi", fil, "-o", d / "mask.json"], d)
        add("fb", "rfi mask", dt, ok)
        for method in ("subband", "fdmt"):
            dt, ok = run([FB, "sp", fil, "-o", d / f"sp_{method}", "--dm", f"0:{a.dm_hi}", "--dedisp", method, "--format", "csv,presto"], d)
            add("fb", f"single-pulse ({method})", dt, ok, read_singlepulse([d / f"sp_{method}.singlepulse"]) if ok else None)
        dt, ok = run([FB, "periodic", fil, "-o", d / "per", "--dm", f"0:{2 * pulsar[0]}", "--fold", "5"], d)
        found = None
        if ok:
            lines = (d / "per.periodic.csv").read_text().splitlines()[1:]
            if lines:
                v = lines[0].split(",")
                found = (float(v[7]), float(v[9]), float(v[2]))
        add("fb", "periodic (fft, fold 5)", dt, ok, found)

    if "presto" not in skip and shutil.which("prepsubband"):
        d = work / "presto"
        d.mkdir(exist_ok=True)
        dt, ok = run(["rfifind", "-time", "2.0", "-o", "bench", fil], d)
        add("presto", "rfifind", dt, ok)
        step = 0.5 if a.dm_hi <= 600 else 1.0
        ndm = int(a.dm_hi / step)
        cmd = ["prepsubband", "-lodm", "0", "-dmstep", str(step), "-numdms", str(min(ndm, 1000)), "-nsub", str(min(128, a.nchan)), "-o", "bench", fil]
        if (d / "bench_rfifind.mask").exists():
            cmd[1:1] = ["-mask", "bench_rfifind.mask"]
        dt, ok = run(cmd, d)
        add("presto", "prepsubband", dt, ok)
        dats = sorted(d.glob("bench_DM*.dat"))
        if dats and shutil.which("single_pulse_search.py"):
            dt, ok = run(["single_pulse_search.py", "-b", "--noplot", *[p.name for p in dats]], d)
            add("presto", "single_pulse_search.py", dt, ok, read_singlepulse(d.glob("bench_DM*.singlepulse")) if ok else None)
        target = min(dats, key=lambda p: abs(float(p.stem.split("DM")[1]) - pulsar[0])) if dats else None
        if target is not None and shutil.which("accelsearch"):
            t0 = time.perf_counter()
            ok = run(["realfft", target.name], d)[1]
            red = target.with_name(target.stem + "_red.fft")
            ok = ok and run(["rednoise", target.with_suffix(".fft").name], d)[1]
            if ok and not red.with_suffix(".inf").exists():
                shutil.copy(target.with_suffix(".inf"), red.with_suffix(".inf"))
            ok = ok and run(["accelsearch", "-zmax", "0", "-numharm", "16", "-flo", "0.5", red.name], d)[1]
            add("presto", "realfft + rednoise + accelsearch (one DM)", time.perf_counter() - t0, ok)

    if "heimdall" not in skip and shutil.which("heimdall"):
        d = work / "heimdall"
        d.mkdir(exist_ok=True)
        dt, ok = run(["heimdall", "-f", fil, "-dm", "0", str(a.dm_hi), "-output_dir", d], d)
        add("heimdall", "single-pulse (GPU)", dt, ok, read_heimdall(d.glob("*.cand")) if ok else None)

    print(f"\nexpected: burst DM {burst[0]} at {burst[1]} s; pulsar DM {pulsar[0]} P {pulsar[1]} s\n")
    print("| tool | step | seconds | result |")
    print("|---|---|---|---|")
    for r in results:
        f = r["found"]
        res = "failed" if not r["ok"] else ("" if f is None else f"DM {f[0]:.1f}, S/N {f[1]:.1f}, {f[2]:.4f}")
        print(f"| {r['tool']} | {r['step']} | {r['seconds']:.2f} | {res} |")
    (work / "results.json").write_text(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
