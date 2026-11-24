import ctypes as C
import json
import os
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parents[1]
_EXT = {"darwin": "dylib", "win32": "dll"}.get(sys.platform, "so")
_NAME = ("" if sys.platform == "win32" else "lib") + f"fbcore.{_EXT}"


class FbError(RuntimeError):
    pass


def _find():
    cands = [os.environ.get("FB_LIB"), _HERE / _NAME, _ROOT / "target" / "release" / _NAME, _ROOT / "target" / "debug" / _NAME]
    for c in cands:
        if c and Path(c).exists():
            return str(c)
    raise FbError(f"{_NAME} not found; run `make` (or `cargo build --release`) in {_ROOT}")


lib = C.CDLL(_find())

_s = C.c_char_p
_p = C.c_void_p
_f32p = C.POINTER(C.c_float)
_f64p = C.POINTER(C.c_double)
_u8p = C.POINTER(C.c_uint8)
_i32p = C.POINTER(C.c_int32)

_SIGS = {
    "fb_version": ([], C.c_char_p),
    "fb_last_error": ([], C.c_char_p),
    "fb_free_string": ([_p], None),
    "fb_open": ([_s, _s], _p),
    "fb_close": ([_p], None),
    "fb_header": ([_p], _p),
    "fb_read": ([_p, C.c_uint64, C.c_uint64, _f32p, C.c_uint64], C.c_int64),
    "fb_write_fil": ([_s, _s, C.c_uint32, _f32p, C.c_uint64], C.c_int64),
    "fb_plan": ([_s, C.c_double, C.c_double, C.c_double, C.c_double, C.c_double], _p),
    "fb_plan_from_list": ([_s, _s], _p),
    "fb_delay_table": ([_s, _i32p, C.c_uint64], C.c_int64),
    "fb_dedisperse": ([_s, _s, _f32p, C.c_uint64, _u8p, C.c_uint32, _f32p, C.c_uint64], C.c_int64),
    "fb_iqrm": ([_f64p, C.c_uint64, C.c_uint32, C.c_double, _u8p], C.c_int64),
    "fb_mask_new": ([_s, C.c_double], _p),
    "fb_rfi_block": ([_s, _f32p, C.c_uint64, C.c_uint64, _s, _s], _p),
    "fb_mask_apply": ([_s, _f32p, C.c_uint64, C.c_uint64, _s, C.c_uint8], C.c_int64),
    "fb_presto_mask_read": ([_s, _s], _p),
    "fb_presto_mask_write": ([_s, C.c_double, _s], C.c_int64),
    "fb_dat_info": ([_s], _p),
    "fb_dat_read": ([_s, _f32p, C.c_uint64], C.c_int64),
    "fb_dat_write": ([_s, _s, C.c_double, C.c_double, _f32p, C.c_uint64], C.c_int64),
    "fb_sp_search": ([_s, _s, _f32p, C.c_uint64, C.c_uint64, C.c_uint64, _s], _p),
    "fb_run_sp": ([_s, _s], _p),
}
for _n, (_a, _r) in _SIGS.items():
    _fn = getattr(lib, _n)
    _fn.argtypes = _a
    _fn.restype = _r


def err():
    return FbError(lib.fb_last_error().decode())


def enc(x):
    if x is None:
        return None
    if isinstance(x, (bytes, str)):
        return x.encode() if isinstance(x, str) else x
    return json.dumps(x).encode()


def take_json(ptr):
    if not ptr:
        raise err()
    try:
        return json.loads(C.string_at(ptr).decode())
    finally:
        lib.fb_free_string(ptr)


def check(v):
    if v < 0:
        raise err()
    return v


def f32(a):
    a = np.ascontiguousarray(a, dtype=np.float32)
    return a, a.ctypes.data_as(_f32p)


def version():
    return lib.fb_version().decode()
