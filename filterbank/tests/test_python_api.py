import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
import filterbank as fb

fails = 0


def check(name, cond):
    global fails
    print(("ok   " if cond else "FAIL ") + name)
    fails += 0 if cond else 1


d = Path(tempfile.mkdtemp())
rng = np.random.default_rng(1)
nchan, nsamp, tsamp = 64, 40000, 256e-6
hdr = fb.Header({
    "source_name": "PYTEST", "telescope": "GMRT", "backend": "FAKE", "ra_deg": 10.0, "dec_deg": -0.25,
    "az_deg": None, "za_deg": None, "start": {"imjd": 60100, "sec": 1234.5}, "tsamp": tsamp, "nsamples": nsamp,
    "freqs": {"mhz": list(1450.0 - 2.0 * np.arange(nchan)), "chan_bw": -2.0}, "npol": 1, "pol_type": "Intensity",
    "nbits": 32, "scaling": {"zero_off": 0.0, "signed": False, "applied": False}, "barycentric": False, "beam": 1,
    "nbeams": 1, "order": "TimePolFreq", "provenance": [],
})
data = rng.normal(20, 2, (nsamp, nchan)).astype(np.float32)
dm0, t0 = 120.0, 4.0
delays = fb.K_DM * dm0 * (hdr.freqs ** -2 - hdr.freqs[0] ** -2)
idx = ((t0 + delays) / tsamp).astype(int)
for c in range(nchan):
    data[idx[c]:idx[c] + 6, c] += 3.0
data[:, 10] += np.abs(rng.normal(0, 30, nsamp))

fb.write_fil(d / "a.fil", hdr, data, nbits=32)
with fb.open(d / "a.fil") as r:
    back = r.read()
    check("fil write/read round trip", np.array_equal(back[:, 0, :], data) and r.header.nchan == nchan)
    check("header epoch kept", abs(r.header.mjd - hdr.mjd) * 86400 < 1e-5 and r.header["dec_deg"] < 0)
    h = r.header

m = fb.rfi.find(h, back, interval=1.0)
check("IQRM masks RFI channel 10", m.channel_fraction()[10] > 0.9 and m.fraction < 0.1)
fb.write_presto_mask(m, h.tsamp, d / "a_rfifind.mask")
m2 = fb.read_presto_mask(d / "a_rfifind.mask", h)
check("PRESTO mask round trip", np.array_equal(m.grid, m2.grid))
cleaned, nzap = fb.rfi.apply(h, back, m)
check("mask applied", nzap > 0 and np.allclose(cleaned[:, 0, 10], cleaned[:, 0, 10].mean()))

p = fb.dedisp.plan(h, 0, 300, step=2.0)
tr = fb.dedisp.dedisperse(h, p, cleaned, killmask=m.channel_fraction() > 0.5)
best = np.unravel_index(np.argmax(tr), tr.shape)
check("dedispersion peak at injected DM", abs(p.dms[best[0]] - dm0) <= 4 and abs(best[1] * tsamp - t0) < 0.01)
tr_sub = fb.dedisp.dedisperse(h, p, cleaned, nsub=8)
check("subband agrees with brute force peak", abs(np.argmax(tr_sub.max(axis=1)) - best[0]) <= 2)

c = fb.singlepulse.search(h, p, tr, threshold=7)
top = max(c, key=lambda x: x["snr"])
check("search finds pulse", abs(top["dm"] - dm0) <= 4 and abs(fb.singlepulse.candidate_time(top, h) - t0) < 0.01)

res = fb.singlepulse.run(d / "a.fil", dm_lo=0, dm_hi=300, dm_step=2.0, threshold=7)
top2 = max(res["candidates"], key=lambda x: x["snr"])
check("pipeline run matches manual chain", abs(top2["dm"] - top["dm"]) <= 2 and len(res["candidates"]) >= 1)

s = tr[best[0]]
fb.write_dat(d / "a_DM120", h, p.dms[best[0]], s)
hh, dmv, fref, s2 = fb.read_dat(d / "a_DM120.dat")
check("PRESTO .dat/.inf round trip", np.array_equal(s, s2) and abs(dmv - p.dms[best[0]]) < 1e-9)
check("backend registry", fb.backend.available()[0] == "cpu")
print(f"{fails} failure(s), library {fb.version()}")
sys.exit(1 if fails else 0)
