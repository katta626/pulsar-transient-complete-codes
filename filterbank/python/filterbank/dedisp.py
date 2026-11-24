import numpy as np

from . import _lib as L
from . import backend as B
from .types import Header, Plan


def plan(header, lo=0.0, hi=1000.0, step=None, width=1e-3, tol=1.25, dms=None):
    h = Header(header)
    if dms is not None:
        return Plan(L.take_json(L.lib.fb_plan_from_list(L.enc(h.core()), L.enc([float(d) for d in dms]))))
    return Plan(L.take_json(L.lib.fb_plan(L.enc(h.core()), lo, hi, step or 0.0, width, tol)))


def delay_table(p):
    nchan = len(p["freqs"]["mhz"])
    out = np.empty((len(p["dms"]), nchan), dtype=np.int32)
    maxd = L.check(L.lib.fb_delay_table(L.enc(dict(p)), out.ctypes.data_as(L._i32p), out.size))
    return out, maxd


def dedisperse(header, p, data, killmask=None, nsub=0, backend="auto"):
    h = Header(header).intensity()
    d = data[:, 0, :] if data.ndim == 3 else data
    d, ptr = L.f32(d)
    be = B.resolve(backend)
    table, maxd = delay_table(p)
    nout = d.shape[0] - maxd
    if nout <= 0:
        raise ValueError(f"need more than {maxd} samples for this plan")
    if be == "cuda":
        from . import cuda as CU
        return CU.dedisperse(d, table, maxd, killmask)
    out = np.empty((len(p["dms"]), nout), dtype=np.float32)
    km = None
    if killmask is not None:
        kmask = np.ascontiguousarray(killmask, dtype=np.uint8)
        km = kmask.ctypes.data_as(L._u8p)
    n = L.check(L.lib.fb_dedisperse(L.enc(h.core()), L.enc(dict(p)), ptr, d.shape[0], km, nsub, out.ctypes.data_as(L._f32p), out.size))
    return out[:, :n]
