import ctypes as C
import os
import sys
from pathlib import Path

import numpy as np

from . import _lib as L
from .types import Header

_EXT = {"darwin": "dylib"}.get(sys.platform, "so")
_ROOT = Path(__file__).resolve().parents[2]
_dada = None


def _load():
    global _dada
    if _dada is not None:
        return _dada
    cands = [os.environ.get("FB_DADA_LIB"), _ROOT / "dada" / "build" / f"libfbdada.{_EXT}"]
    for c in cands:
        if c and Path(c).exists():
            lib = C.CDLL(str(c))
            break
    else:
        raise L.FbError(f"libfbdada.{_EXT} not found; run `make dada` in {_ROOT}")
    sig = {
        "fbdada_create": ([C.c_int, C.c_uint64, C.c_uint64, C.c_uint64, C.c_uint64, C.c_uint], C.c_int),
        "fbdada_destroy": ([C.c_int], C.c_int),
        "fbdada_open_reader": ([C.c_int], C.c_void_p),
        "fbdada_open_writer": ([C.c_int], C.c_void_p),
        "fbdada_next_header": ([C.c_void_p], C.c_char_p),
        "fbdada_read": ([C.c_void_p, C.c_char_p, C.c_uint64], C.c_int64),
        "fbdada_write_header": ([C.c_void_p, C.c_char_p], C.c_int),
        "fbdada_write": ([C.c_void_p, C.c_char_p, C.c_uint64], C.c_int64),
        "fbdada_bufsz": ([C.c_void_p], C.c_uint64),
        "fbdada_close": ([C.c_void_p], None),
    }
    for n, (a, r) in sig.items():
        f = getattr(lib, n)
        f.argtypes = a
        f.restype = r
    _dada = lib
    return lib


def available():
    try:
        _load()
        return True
    except L.FbError:
        return False


def _key(key):
    return int(key, 16) if isinstance(key, str) else int(key)


def _sysv_limit():
    if sys.platform != "darwin":
        return None
    import subprocess

    try:
        v = subprocess.run(["sysctl", "-n", "kern.sysv.shmmax", "kern.sysv.shmall", "kern.sysv.shmseg"], capture_output=True, text=True).stdout.split()
        return min(int(v[0]), int(v[1]) * 4096), int(v[2])
    except Exception:
        return None


def create(key, nbufs=8, bufsz=1 << 22, hdr_nbufs=8, hdr_bufsz=4096, nreaders=1):
    lim = _sysv_limit()
    if lim is not None:
        limit, nseg = lim
        total = nbufs * bufsz + hdr_nbufs * hdr_bufsz + 4096
        if total > limit:
            raise L.FbError(f"ring needs {total} bytes of System V shared memory but this Mac allows {limit}; raise kern.sysv.shmmax/shmall or use a smaller ring")
        if nbufs + hdr_nbufs + 2 > nseg:
            raise L.FbError(f"ring needs {nbufs + hdr_nbufs + 2} shared-memory segments per process but this Mac allows {nseg} (kern.sysv.shmseg)")
    if _load().fbdada_create(_key(key), nbufs, bufsz, hdr_nbufs, hdr_bufsz, nreaders) != 0:
        raise L.FbError(f"could not create ring {key}")


def destroy(key):
    return _load().fbdada_destroy(_key(key)) == 0


class Ring:
    def __init__(self, key, mode="r"):
        lib = _load()
        self.lib = lib
        self.mode = mode
        opener = lib.fbdada_open_reader if mode == "r" else lib.fbdada_open_writer
        self.h = opener(_key(key))
        if not self.h:
            raise L.FbError(f"could not attach to ring {key} ({mode})")

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()

    def close(self):
        if self.h:
            self.lib.fbdada_close(self.h)
            self.h = None

    def header_text(self):
        t = self.lib.fbdada_next_header(self.h)
        if t is None:
            raise L.FbError("no header in ring")
        return t.decode(errors="replace")

    def read(self, nbytes):
        buf = C.create_string_buffer(nbytes)
        n = self.lib.fbdada_read(self.h, buf, nbytes)
        if n < 0:
            raise L.FbError("ring read failed")
        return buf.raw[:n]

    def write_header(self, text):
        if self.lib.fbdada_write_header(self.h, text.encode()) < 0:
            raise L.FbError("ring header write failed")

    def write(self, data):
        data = bytes(data)
        if self.lib.fbdada_write(self.h, data, len(data)) != len(data):
            raise L.FbError("ring write failed")


def layout(text, sum_pols=True):
    info = L.take_json(L.lib.fb_dada_layout(text.encode(), int(sum_pols)))
    return Header(info["header"]), info


def decode(text, raw, sum_pols=True, info=None):
    if info is None:
        _, info = layout(text, sum_pols)
    h = info["header"]
    nchan = len(h["freqs"]["mhz"])
    npol = h["npol"]
    nsamp = len(raw) // info["bytes_per_sample"]
    out = np.empty(nsamp * npol * nchan, dtype=np.float32)
    src = np.frombuffer(raw, dtype=np.uint8)
    n = L.lib.fb_dada_decode(text.encode(), int(sum_pols), src.ctypes.data_as(L._u8p), len(raw), out.ctypes.data_as(L._f32p), out.size)
    if n < 0:
        raise L.err()
    return out.reshape(nsamp, npol, nchan)


def follow(key, gulp=65536, overlap=0, sum_pols=True, on_header=None):
    with Ring(key, "r") as r:
        text = r.header_text()
        header, info = layout(text, sum_pols)
        if on_header is not None:
            ov = on_header(header)
            overlap = overlap if ov is None else int(ov)
        bps = info["bytes_per_sample"]
        pos = 0
        carry = b""
        need = (gulp + overlap) * bps
        while True:
            chunk = r.read(need - len(carry))
            eod = len(chunk) < need - len(carry)
            buf = carry + chunk
            if len(buf) >= need or (eod and len(buf) > overlap * bps):
                yield pos, decode(text, buf, sum_pols, info)
                keep = overlap * bps if len(buf) >= need else 0
                carry = buf[len(buf) - keep:] if keep else b""
                pos += (len(buf) - keep) // bps
            if eod:
                return


def play(path, key, chunk=1 << 20):
    path = Path(path)
    with open(path, "rb") as f:
        head = f.read(65536)
    text = head.split(b"\0", 1)[0].decode(errors="replace")
    hdr_size = 4096
    for line in text.splitlines():
        p = line.split()
        if len(p) >= 2 and p[0].upper() == "HDR_SIZE":
            hdr_size = int(p[1])
    text = head[:hdr_size].split(b"\0", 1)[0].decode(errors="replace")
    with Ring(key, "w") as w, open(path, "rb") as f:
        w.write_header(text)
        f.seek(hdr_size)
        while True:
            b = f.read(chunk)
            if not b:
                break
            w.write(b)
