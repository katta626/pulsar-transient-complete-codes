import numpy as np

from . import dedisp
from .io import Reader
from .types import K_DM, Header


def _resize(a, shape):
    out = a
    for axis, n in enumerate(shape):
        m = out.shape[axis]
        if m == n:
            continue
        if m > n and m % n == 0:
            new = list(out.shape)
            new[axis] = n
            new.insert(axis + 1, m // n)
            out = out.reshape(new).mean(axis=axis + 1)
        else:
            x_old = np.linspace(0, 1, m)
            x_new = np.linspace(0, 1, n)
            out = np.apply_along_axis(lambda v: np.interp(x_new, x_old, v), axis, out)
    return out


def _norm(a):
    a = a - np.median(a)
    s = 1.4826 * np.median(np.abs(a)) or a.std() or 1.0
    return (a / s).astype(np.float32)


def cutout(reader, cand, nt=256, nf=256, ndm=256, dm_span=None):
    h = reader.header
    tsamp = h.tsamp
    width = max(1, int(cand["width_samples"]))
    dm = float(cand["dm"])
    freqs = h.freqs
    sweep = K_DM * dm * (freqs.min() ** -2 - freqs.max() ** -2)
    dec = max(1, width // 2)
    half = nt * dec // 2
    centre = int(cand["sample"])
    start = max(0, centre - half)
    extra = int(np.ceil(K_DM * max(2 * dm, dm + 50) * (freqs.min() ** -2 - freqs.max() ** -2) / tsamp)) + 2
    block = reader.read(start, 2 * half + extra)[:, 0, :]
    if block.shape[0] < 2 * half + 1:
        block = np.pad(block, ((0, 2 * half + extra - block.shape[0]), (0, 0)), mode="edge")
    delays = np.round(K_DM * dm * (freqs ** -2 - freqs.max() ** -2) / tsamp).astype(int)
    ft = np.stack([block[delays[c]:delays[c] + 2 * half, c] for c in range(len(freqs))], axis=1)
    ft = ft[: (ft.shape[0] // dec) * dec].reshape(-1, dec, ft.shape[1]).mean(axis=1)
    ft = _resize(ft, (nt, nf)).T
    span = dm_span if dm_span is not None else max(2 * dm, 20.0)
    dms = np.linspace(max(0.0, dm - span / 2), dm + span / 2, ndm)
    plan = dedisp.plan(h, dms=dms)
    table, maxd = dedisp.delay_table(plan)
    need = 2 * half + maxd
    blk = block if block.shape[0] >= need else np.pad(block, ((0, need - block.shape[0]), (0, 0)), mode="edge")
    dmt = dedisp.dedisperse(h, plan, blk[:need], backend="cpu")[:, : 2 * half]
    dmt = dmt[:, : (dmt.shape[1] // dec) * dec].reshape(ndm, -1, dec).mean(axis=2)
    dmt = _resize(dmt, (ndm, nt))
    return {"freq_time": _norm(ft), "dm_time": _norm(dmt), "dms": dms, "sweep_s": sweep}


def features(c, nsub=16):
    ft = c["freq_time"]
    dmt = c["dm_time"]
    nf, nt = ft.shape
    prof = ft.mean(axis=0)
    centre = nt // 2
    lo, hi = max(0, centre - nt // 8), min(nt, centre + nt // 8)
    peak_t = lo + int(np.argmax(prof[lo:hi]))
    above = prof > 0.5 * prof[peak_t]
    a = peak_t
    while a > 0 and above[a - 1]:
        a -= 1
    b = peak_t + 1
    while b < nt and above[b]:
        b += 1
    on = slice(a, b)
    off = np.r_[0:max(1, a - 4), min(nt, b + 4):nt]
    sub = ft[: (nf // nsub) * nsub].reshape(nsub, -1, nt).mean(axis=1)
    sub_on = sub[:, on].sum(axis=1)
    sub_noise = sub[:, off].std(axis=1) * np.sqrt(b - a) + 1e-6
    broadband = float(np.mean(sub_on / sub_noise > 1.5))
    curve = dmt[:, on].mean(axis=1)
    pk = int(np.argmax(curve))
    zero_dm_ratio = float(max(curve[0], 0) / (curve[pk] + 1e-6))
    dm_peak_width = float(np.mean(curve > 0.5 * curve[pk]))
    offset = abs(pk - len(curve) / 2) / len(curve)
    return {"broadband": broadband, "zero_dm_ratio": zero_dm_ratio, "dm_peak_width": dm_peak_width, "dm_peak_offset": float(offset), "on_width": int(b - a)}


def heuristic_score(f):
    s = 1.0
    s *= np.clip((f["broadband"] - 0.15) / 0.45, 0, 1)
    s *= np.clip(1.0 - f["zero_dm_ratio"], 0, 1)
    s *= np.clip(1.0 - f["dm_peak_offset"] * 4, 0, 1)
    s *= np.clip(1.2 - f["dm_peak_width"], 0, 1)
    return float(s)


class OnnxModel:
    def __init__(self, path):
        try:
            import onnxruntime as ort
        except ImportError as e:
            raise RuntimeError("onnxruntime is not installed; use heuristic scoring or install it") from e
        self.sess = ort.InferenceSession(str(path))
        self.inputs = [i.name for i in self.sess.get_inputs()]

    def __call__(self, c):
        ft = c["freq_time"][None, :, :, None].transpose(0, 2, 1, 3)
        dmt = c["dm_time"][None, :, :, None]
        feeds = dict(zip(self.inputs, [ft, dmt][: len(self.inputs)]))
        out = self.sess.run(None, feeds)[0]
        return float(np.ravel(out)[-1])


def classify(paths, candidates, model=None, min_score=0.5, save_dir=None, **kw):
    r = paths if isinstance(paths, Reader) else Reader(paths)
    m = OnnxModel(model) if model else None
    out = []
    for i, cand in enumerate(candidates):
        c = cutout(r, cand, **kw)
        f = features(c)
        score = m(c) if m else heuristic_score(f)
        rec = dict(cand)
        rec.update(f)
        rec["probability"] = score
        rec["label"] = "astro" if score >= min_score else "rfi"
        out.append(rec)
        if save_dir is not None:
            from pathlib import Path
            p = Path(save_dir)
            p.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(p / f"cand_{i:05d}_dm{cand['dm']:.2f}_s{cand['sample']}.npz", **c, **{k: v for k, v in rec.items() if isinstance(v, (int, float, str))})
    return out
