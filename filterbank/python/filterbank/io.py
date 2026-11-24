import numpy as np

from . import _lib as L
from .types import Header, Mask


class Reader:
    def __init__(self, paths, *, zero_off=None, scales=None, weights=None, sum_pols=True, raw=False):
        if isinstance(paths, (str, bytes)) or hasattr(paths, "__fspath__"):
            paths = [paths]
        opts = {"sum_pols": sum_pols, "raw": raw}
        if zero_off:
            opts["zero_off"] = zero_off
        if scales:
            opts["scales"] = scales
        if weights is not None:
            opts["weights"] = weights
        self.paths = [str(p) for p in paths]
        self._h = L.lib.fb_open(L.enc(self.paths), L.enc(opts))
        if not self._h:
            raise L.err()
        self.header = Header(L.take_json(L.lib.fb_header(self._h)))

    def read(self, start=0, nsamp=None):
        h = self.header
        if nsamp is None:
            nsamp = h.nsamples - start
        nsamp = max(0, min(nsamp, h.nsamples - start))
        out = np.empty((nsamp, h.npol, h.nchan), dtype=np.float32)
        n = L.check(L.lib.fb_read(self._h, start, nsamp, out.ctypes.data_as(L._f32p), out.size))
        return out[:n]

    def blocks(self, gulp=65536, overlap=0):
        s = 0
        while s < self.header.nsamples:
            b = self.read(s, gulp + overlap)
            if len(b) == 0:
                break
            yield s, b
            s += gulp

    def close(self):
        if self._h:
            L.lib.fb_close(self._h)
            self._h = None

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass


def open(paths, **kw):
    return Reader(paths, **kw)


def write_fil(path, header, data, nbits=32):
    h = Header(header)
    data, ptr = L.f32(data)
    return L.check(L.lib.fb_write_fil(str(path).encode(), L.enc(h.core()), nbits, ptr, data.size))


def read_dat(path):
    info = L.take_json(L.lib.fb_dat_info(str(path).encode()))
    out = np.empty(info["n"], dtype=np.float32)
    L.check(L.lib.fb_dat_read(str(path).encode(), out.ctypes.data_as(L._f32p), out.size))
    return Header(info["header"]), info["dm"], info["f_ref"], out


def write_dat(root, header, dm, data, f_ref=None):
    h = Header(header)
    f_ref = float(max(h.freqs)) if f_ref is None else f_ref
    data, ptr = L.f32(data)
    return L.check(L.lib.fb_dat_write(str(root).encode(), L.enc(h.core()), dm, f_ref, ptr, data.size))


def read_presto_mask(path, header=None):
    return Mask(L.take_json(L.lib.fb_presto_mask_read(str(path).encode(), L.enc(Header(header).core()) if header else None)))


def write_presto_mask(mask, tsamp, path):
    return L.check(L.lib.fb_presto_mask_write(L.enc(dict(mask)), tsamp, str(path).encode()))
