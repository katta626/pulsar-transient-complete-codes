import ctypes as C

import numpy as np

from . import backend as B

_f32p = C.POINTER(C.c_float)
_i32p = C.POINTER(C.c_int32)
_u8p = C.POINTER(C.c_uint8)


def _lib():
    lib = B.cuda()
    if lib is None:
        raise RuntimeError("CUDA backend not available")
    if not getattr(lib, "_fb_typed", False):
        lib.fbcuda_last_error.restype = C.c_char_p
        lib.fbcuda_dedisperse.argtypes = [_f32p, C.c_int, C.c_int, _i32p, C.c_int, _u8p, C.c_int, _f32p]
        lib.fbcuda_zero_dm.argtypes = [_f32p, C.c_int, C.c_int, C.c_int]
        lib.fbcuda_boxcar_snr.argtypes = [_f32p, C.c_int, C.c_int, _i32p, C.c_int, _f32p, _i32p]
        lib.fbcuda_set_device.argtypes = [C.c_int]
        lib._fb_typed = True
    return lib


def _check(lib, rc):
    if rc != 0:
        raise RuntimeError(lib.fbcuda_last_error().decode())


def set_device(i):
    lib = _lib()
    _check(lib, lib.fbcuda_set_device(i))


def dedisperse(data, delay_table, max_delay, killmask=None):
    lib = _lib()
    data = np.ascontiguousarray(data, dtype=np.float32)
    table = np.ascontiguousarray(delay_table, dtype=np.int32)
    nsamp, nchan = data.shape
    ndm = table.shape[0]
    nout = nsamp - max_delay
    out = np.empty((ndm, nout), dtype=np.float32)
    km = None
    if killmask is not None:
        k = np.ascontiguousarray(killmask, dtype=np.uint8)
        km = k.ctypes.data_as(_u8p)
    _check(lib, lib.fbcuda_dedisperse(data.ctypes.data_as(_f32p), nsamp, nchan, table.ctypes.data_as(_i32p), ndm, km, nout, out.ctypes.data_as(_f32p)))
    return out


def zero_dm(data, keep_channel_mean=True):
    lib = _lib()
    d = np.ascontiguousarray(data, dtype=np.float32).copy()
    _check(lib, lib.fbcuda_zero_dm(d.ctypes.data_as(_f32p), d.shape[0], d.shape[1], int(keep_channel_mean)))
    return d


def boxcar_snr(series, widths):
    lib = _lib()
    s = np.ascontiguousarray(series, dtype=np.float32)
    w = np.ascontiguousarray(widths, dtype=np.int32)
    ndm, nout = s.shape
    snr = np.empty_like(s)
    bw = np.empty(s.shape, dtype=np.int32)
    _check(lib, lib.fbcuda_boxcar_snr(s.ctypes.data_as(_f32p), ndm, nout, w.ctypes.data_as(_i32p), len(w), snr.ctypes.data_as(_f32p), bw.ctypes.data_as(_i32p)))
    return snr, bw


def _typed_extra(lib):
    if not getattr(lib, "_fb_typed2", False):
        lib.fbcuda_fdmt.argtypes = [_f32p, C.c_int, C.c_int, _i32p, C.c_int, C.c_int, C.c_int, _i32p, _i32p, _i32p, _i32p, C.c_int, C.c_int, _i32p, _i32p, C.c_int, _f32p]
        lib.fbcuda_fft_search.argtypes = [_f32p, C.c_int, C.c_int, C.c_int, _i32p, _i32p, _i32p, _f32p, C.c_int, _i32p, _i32p, _i32p, _f32p]
        lib.fbcuda_fft_search.restype = C.c_int
        lib._fb_typed2 = True
    return lib


def _i32(a):
    a = np.ascontiguousarray(a, dtype=np.int32)
    return a, a.ctypes.data_as(_i32p)


def fdmt(header, p, data, killmask=None):
    from . import _lib as L
    from .types import Header

    lib = _typed_extra(_lib())
    h = Header(header).intensity()
    d = np.ascontiguousarray(data[:, 0, :] if data.ndim == 3 else data, dtype=np.float32)
    nsamp, nchan = d.shape
    km = None
    if killmask is not None:
        kb = np.ascontiguousarray(killmask, dtype=np.uint8)
        km = kb.ctypes.data_as(L._u8p)
    s = L.take_json(L.lib.fb_fdmt_setup(L.enc(h.core()), L.enc(dict(p)), nsamp, km))
    levels = s["levels"]
    nrows, nrows_p = _i32([len(l["lo"]) for l in levels])
    lo, lo_p = _i32(np.concatenate([l["lo"] for l in levels]) if levels else [0])
    hi, hi_p = _i32(np.concatenate([l["hi"] for l in levels]) if levels else [0])
    sh, sh_p = _i32(np.concatenate([l["shift"] for l in levels]) if levels else [0])
    order, order_p = _i32(s["order"])
    rows, rows_p = _i32(s["out_rows"])
    oshift, oshift_p = _i32(s["out_shift"])
    nout = s["nout"]
    out = np.empty((len(rows), nout), dtype=np.float32)
    _check(lib, lib.fbcuda_fdmt(d.ctypes.data_as(_f32p), nsamp, nchan, order_p, s["nf"], s["dt0"] + 1, len(levels), nrows_p, lo_p, hi_p, sh_p, s["max_rows"], len(rows), rows_p, oshift_p, nout, out.ctypes.data_as(_f32p)))
    return out


def _power_threshold(sigma, numharm, ntrials):
    from . import _lib as L

    p = float(numharm)
    while L.lib.fb_sigma_from_power(p, numharm, ntrials) < sigma and p < 1e6:
        p *= 1.1
    return p


def fft_search(series, tsamp, f_min=0.1, f_max=2000.0, max_harm=16, sigma=6.0, cap=1 << 20):
    from . import _lib as L

    lib = _typed_extra(_lib())
    x = np.ascontiguousarray(np.atleast_2d(series), dtype=np.float32)
    nrows, n = x.shape
    if n % 2:
        x = np.ascontiguousarray(x[:, :-1])
        n -= 1
    nbins = n // 2 + 1
    t_obs = n * tsamp
    ntrials = nbins * 2.0
    harms = [1]
    while harms[-1] * 2 <= max_harm:
        harms.append(harms[-1] * 2)
    lo = [max(1, int(f_min * t_obs * 2 * h)) for h in harms]
    hi = [min(int(f_max * t_obs * 2 * h), 2 * nbins - 1) for h in harms]
    thr = [_power_threshold(sigma, h, ntrials) for h in harms]
    harms_a, harms_p = _i32(harms)
    lo_a, lo_p = _i32(lo)
    hi_a, hi_p = _i32(hi)
    thr_a = np.ascontiguousarray(thr, dtype=np.float32)
    c_row = np.empty(cap, np.int32)
    c_stage = np.empty(cap, np.int32)
    c_bin = np.empty(cap, np.int32)
    c_pow = np.empty(cap, np.float32)
    count = lib.fbcuda_fft_search(x.ctypes.data_as(_f32p), nrows, n, len(harms), harms_p, lo_p, hi_p, thr_a.ctypes.data_as(_f32p), cap,
                                  c_row.ctypes.data_as(_i32p), c_stage.ctypes.data_as(_i32p), c_bin.ctypes.data_as(_i32p), c_pow.ctypes.data_as(_f32p))
    if count < 0:
        _check(lib, count)
    m = min(count, cap)
    order = np.lexsort((c_bin[:m], c_stage[:m], c_row[:m]))
    out = [[] for _ in range(nrows)]
    i = 0
    while i < m:
        j = i
        best = order[i]
        while j + 1 < m and c_row[order[j + 1]] == c_row[order[i]] and c_stage[order[j + 1]] == c_stage[order[i]] and c_bin[order[j + 1]] == c_bin[order[j]] + 1:
            j += 1
            if c_pow[order[j]] > c_pow[best]:
                best = order[j]
        h = harms[c_stage[best]]
        pw = float(c_pow[best])
        out[c_row[best]].append({
            "f_hz": c_bin[best] / 2.0 / t_obs / h, "numharm": h, "power": pw, "snr": float(np.sqrt(2 * pw)),
            "sigma": L.lib.fb_sigma_from_power(pw, h, ntrials), "method": "fft (cuda)",
        })
        i = j + 1
    for row in out:
        row.sort(key=lambda c: -c["sigma"])
    return out
