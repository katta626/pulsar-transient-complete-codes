import numpy as np

from . import _lib as L
from .types import Header, Mask


def iqrm(x, radius=None, threshold=3.0):
    x = np.ascontiguousarray(x, dtype=np.float64)
    radius = radius or max(1, len(x) // 10)
    out = np.zeros(len(x), dtype=np.uint8)
    L.check(L.lib.fb_iqrm(x.ctypes.data_as(L._f64p), len(x), radius, threshold, out.ctypes.data_as(L._u8p)))
    return out.astype(bool)


def new_mask(header, interval=1.0):
    return Mask(L.take_json(L.lib.fb_mask_new(L.enc(Header(header).core()), interval)))


def find(header, data, start_sample=0, mask=None, interval=1.0, **cfg):
    h = Header(header)
    if data.ndim == 3:
        data = data[:, 0, :]
    mask = mask or new_mask(h, interval)
    arr = np.ascontiguousarray(data, dtype=np.float32)
    m = L.take_json(L.lib.fb_rfi_block(L.enc(h.intensity().core()), arr.ctypes.data_as(L._f32p), arr.shape[0], start_sample, L.enc(dict(mask)), L.enc(cfg)))
    if cfg.get("zero_dm"):
        data[...] = arr
    return Mask(m)


def apply(header, data, mask, start_sample=0, fill_mean=True):
    h = Header(header)
    squeeze = data.ndim == 2
    arr = np.ascontiguousarray(data if not squeeze else data[:, None, :], dtype=np.float32)
    n = L.check(L.lib.fb_mask_apply(L.enc(h.core()), arr.ctypes.data_as(L._f32p), arr.shape[0], start_sample, L.enc(dict(mask)), int(fill_mean)))
    return (arr[:, 0, :] if squeeze else arr), n


def zero_dm(data, keep_channel_mean=True):
    d = data if data.ndim == 2 else data[:, 0, :]
    m = d.mean(axis=1, keepdims=True)
    out = d - m
    if keep_channel_mean:
        out += d.mean(axis=0, keepdims=True)
    return out
