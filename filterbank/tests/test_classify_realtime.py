import struct
import sys
import tempfile
import threading
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "tests"))
import filterbank as fb
from test_sp_cli import write_fil

fails = 0


def check(name, cond):
    global fails
    print(("ok   " if cond else "FAIL ") + name)
    fails += 0 if cond else 1


d = Path(tempfile.mkdtemp())
rng = np.random.default_rng(8)
nchan, nsamp, tsamp, fch1, foff = 128, 120_000, 256e-6, 1460.0, -1.5
freqs = fch1 + foff * np.arange(nchan)
data = rng.normal(64, 5, size=(nsamp, nchan))
t = np.arange(nsamp) * tsamp
dm, t0 = 180.0, 12.0
for c in range(nchan):
    data[:, c] += 3.5 * np.exp(-0.5 * ((t - t0 - fb.K_DM * dm * (freqs[c] ** -2 - freqs[0] ** -2)) / 1.5e-3) ** 2)
data[int(20.0 / tsamp):int(20.0 / tsamp) + 8, :] += 15.0
write_fil(d / "c.fil", data, fch1, foff, tsamp)

res = fb.singlepulse.run(d / "c.fil", dm_lo=0, dm_hi=400, threshold=8, max_width=64)
cands = sorted(res["candidates"], key=lambda c: -c["snr"])[:10]
out = fb.classify.classify(d / "c.fil", cands, save_dir=d / "cuts")
astro = [c for c in out if c["label"] == "astro"]
best_astro = max(astro, key=lambda c: c["snr"]) if astro else None
check("classifier keeps the dispersed burst", best_astro is not None and abs(best_astro["dm"] - dm) < 15)
hdr = fb.open(d / "c.fil").header
zero_dm = [c for c in out if 19.9 < fb.singlepulse.candidate_time(c, hdr) < 20.1]
check("classifier rejects the zero-DM broadband spike", all(c["label"] == "rfi" for c in zero_dm) and len(zero_dm) > 0)
check("cut-outs saved", len(list((d / "cuts").glob("*.npz"))) == len(out))
z = np.load(next((d / "cuts").glob("*.npz")))
check("cut-out shapes 256x256", z["freq_time"].shape == (256, 256) and z["dm_time"].shape == (256, 256))

live = d / "live.fil"
hdr_bytes = (d / "c.fil").read_bytes()
hlen = len(hdr_bytes) - nsamp * nchan
body = hdr_bytes[hlen:]


def writer():
    with open(live, "wb") as f:
        f.write(hdr_bytes[:hlen])
        f.flush()
        step = nchan * 20000
        for i in range(0, len(body), step):
            f.write(body[i:i + step])
            f.flush()
            time.sleep(0.2)


th = threading.Thread(target=writer)
th.start()
time.sleep(0.3)
found = fb.realtime.run_live(live, gulp=16384, poll=0.1, idle_timeout=1.5, dm_lo=0, dm_hi=400, threshold=8, max_width=64, backend="cpu")
th.join()
h = fb.open(live).header
hit = [c for c in found if abs(fb.singlepulse.candidate_time(c, h) - t0) < 0.02 and abs(c["dm"] - dm) < 15]
check(f"live follow finds the burst while the file grows ({len(found)} candidates)", len(hit) >= 1)
print(f"{fails} failure(s)")
sys.exit(1 if fails else 0)
