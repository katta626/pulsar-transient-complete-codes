import os
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

try:
    from astropy.io import fits
except ImportError:
    print("astropy not installed; skipping")
    sys.exit(0)

ROOT = Path(__file__).resolve().parents[1]
FB = ROOT / "target" / "release" / "fb"


def read_fil(path):
    buf = Path(path).read_bytes()
    pos = 0

    def rstr():
        nonlocal pos
        n = struct.unpack_from("<i", buf, pos)[0]
        pos += 4
        s = buf[pos:pos + n].decode()
        pos += n
        return s

    assert rstr() == "HEADER_START"
    h = {}
    ints = {"telescope_id", "machine_id", "data_type", "barycentric", "nbits", "nchans", "nifs", "nbeams", "ibeam"}
    while True:
        k = rstr()
        if k == "HEADER_END":
            break
        if k in ("rawdatafile", "source_name"):
            h[k] = rstr()
        elif k in ints:
            h[k] = struct.unpack_from("<i", buf, pos)[0]
            pos += 4
        elif k == "signed":
            h[k] = buf[pos]
            pos += 1
        else:
            h[k] = struct.unpack_from("<d", buf, pos)[0]
            pos += 8
    raw = np.frombuffer(buf[pos:], dtype=np.uint8)
    nb = h["nbits"]
    if nb == 32:
        data = np.frombuffer(buf[pos:], dtype="<f4").astype(np.float64)
    elif nb == 8:
        data = raw.astype(np.float64)
    else:
        per = 8 // nb
        shifts = np.arange(per) * nb
        data = ((raw[:, None] >> shifts[None, :]) & ((1 << nb) - 1)).reshape(-1).astype(np.float64)
    width = h["nchans"] * h.get("nifs", 1)
    return h, data[: (data.size // width) * width].reshape(-1, width)


def unpack_msb(raw, nbits):
    if nbits == 8:
        return raw.astype(np.float64)
    per = 8 // nbits
    shifts = 8 - nbits * (np.arange(per) + 1)
    return ((raw[..., None] >> shifts) & ((1 << nbits) - 1)).reshape(*raw.shape[:-1], -1).astype(np.float64)


def make_psrfits(path, *, nbits=8, nchan=16, npol=1, nsblk=32, nrow=3, freqs=None, scl=None, offs=None, wts=None,
                 zero_off=None, stt_smjd=3600, stt_offs=0.25, offs_sub=None, dec="-00:30:00.0", seed=0, pol_type="AA+BB"):
    rng = np.random.default_rng(seed)
    maxv = (1 << nbits) - 1
    vals = rng.integers(0, maxv + 1, size=(nrow, nsblk, npol, nchan), dtype=np.uint8)
    per = 8 // nbits
    shifts = 8 - nbits * (np.arange(per) + 1)
    grouped = vals.reshape(nrow, nsblk, npol, nchan // per, per).astype(np.uint16)
    packed = np.bitwise_or.reduce(grouped << shifts, axis=-1).astype(np.uint8)
    nbytes = nsblk * npol * nchan * nbits // 8
    p = fits.PrimaryHDU()
    hh = p.header
    hh["OBS_MODE"] = "SEARCH"
    hh["SRC_NAME"] = "J1234-0030"
    hh["TELESCOP"] = "GMRT"
    hh["BACKEND"] = "GWB"
    hh["RA"] = "12:34:56.7"
    hh["DEC"] = dec
    hh["STT_IMJD"] = 60000
    hh["STT_SMJD"] = stt_smjd
    hh["STT_OFFS"] = stt_offs
    hh["OBSFREQ"] = 1400.0
    hh["OBSBW"] = -200.0
    if freqs is None:
        freqs = 1500.0 - 200.0 / nchan * np.arange(nchan)
    tbin = 1e-4
    tsub = nsblk * tbin
    if offs_sub is None:
        offs_sub = (np.arange(nrow) + 0.5) * tsub
    cols = [
        fits.Column("TSUBINT", "D", array=np.full(nrow, tsub)),
        fits.Column("OFFS_SUB", "D", array=np.asarray(offs_sub, float)),
        fits.Column("TEL_AZ", "D", array=np.full(nrow, 12.5)),
        fits.Column("TEL_ZEN", "D", array=np.full(nrow, 30.0)),
        fits.Column("DAT_FREQ", f"{nchan}D", array=np.tile(freqs, (nrow, 1))),
        fits.Column("DAT_WTS", f"{nchan}E", array=np.tile(np.ones(nchan) if wts is None else wts, (nrow, 1))),
        fits.Column("DAT_OFFS", f"{nchan * npol}E", array=np.tile(np.zeros(nchan * npol) if offs is None else offs, (nrow, 1))),
        fits.Column("DAT_SCL", f"{nchan * npol}E", array=np.tile(np.ones(nchan * npol) if scl is None else scl, (nrow, 1))),
        fits.Column("DATA", f"{nbytes}B", array=packed.reshape(nrow, nbytes)),
    ]
    t = fits.BinTableHDU.from_columns(cols, name="SUBINT")
    for k, v in dict(NCHAN=nchan, NPOL=npol, NBITS=nbits, NSBLK=nsblk, TBIN=tbin, CHAN_BW=float(freqs[1] - freqs[0]),
                     POL_TYPE=pol_type, NSUBOFFS=0).items():
        t.header[k] = v
    if zero_off is not None:
        t.header["ZERO_OFF"] = zero_off
    fits.HDUList([p, t]).writeto(path, overwrite=True)
    return vals.reshape(nrow * nsblk, npol, nchan).astype(np.float64)


def fb(*args):
    r = subprocess.run([str(FB), *map(str, args)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r


def test_all():
    d = Path(tempfile.mkdtemp())
    failures = 0

    def check(name, cond):
        nonlocal failures
        print(("ok   " if cond else "FAIL ") + name)
        failures += 0 if cond else 1

    for nbits in (2, 4, 8):
        vals = make_psrfits(d / f"b{nbits}.fits", nbits=nbits)
        fb("convert", d / f"b{nbits}.fits", "-o", d / f"b{nbits}.fil", "--overwrite")
        h, data = read_fil(d / f"b{nbits}.fil")
        check(f"{nbits}-bit lossless, same nbits", h["nbits"] == nbits and np.array_equal(data, vals[:, 0, :]))

    a = make_psrfits(d / "c1.fits", seed=1, stt_smjd=3600, stt_offs=0.0)
    b = make_psrfits(d / "c2.fits", seed=2, stt_smjd=3600, stt_offs=3 * 32 * 1e-4 + 2 * 32 * 1e-4)
    fb("convert", d / "c2.fits", d / "c1.fits", "-o", d / "c.fil", "--overwrite")
    h, data = read_fil(d / "c.fil")
    expect = np.concatenate([a[:, 0, :], np.zeros((64, 16)), b[:, 0, :]])
    check("files sorted by time, 64-sample gap filled", data.shape == expect.shape and np.array_equal(data, expect))

    rng = np.random.default_rng(5)
    scl = rng.uniform(0.5, 2.0, 16)
    offs = rng.uniform(-3, 3, 16)
    wts = np.ones(16)
    wts[4] = 0.0
    v = make_psrfits(d / "s.fits", scl=scl, offs=offs, wts=wts, zero_off=127.5)
    fb("convert", d / "s.fits", "-o", d / "s.fil", "--overwrite")
    h, data = read_fil(d / "s.fil")
    expect = ((v[:, 0, :] - 127.5) * scl + offs) * wts
    check("scales/offsets/weights/ZERO_OFF applied, float32 out", h["nbits"] == 32 and np.allclose(data, expect, atol=1e-4))

    v = make_psrfits(d / "p.fits", npol=2, pol_type="AABB")
    fb("convert", d / "p.fits", "-o", d / "p.fil", "--overwrite")
    h, data = read_fil(d / "p.fil")
    check("AABB summed to intensity, nifs=1", h["nifs"] == 1 and np.array_equal(data, v[:, 0, :] + v[:, 1, :]))

    asc = 1300.0 + 200.0 / 16 * np.arange(16)
    v = make_psrfits(d / "a.fits", freqs=asc)
    fb("convert", d / "a.fits", "-o", d / "a.fil", "--overwrite")
    h, data = read_fil(d / "a.fil")
    check("ascending band flipped to descending", h["foff"] < 0 and abs(h["fch1"] - asc[-1]) < 1e-9 and np.array_equal(data, v[:, 0, ::-1]))

    check("negative declination kept", h["src_dej"] < 0 and abs(h["src_dej"] + 3000.0) < 1e-6)
    check("start time precise", abs(h["tstart"] - (60000 + (3600 + 0.25) / 86400)) < 1e-11)

    print(f"{failures} failure(s)")
    return failures


if __name__ == "__main__":
    sys.exit(1 if test_all() else 0)
