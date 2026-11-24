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
