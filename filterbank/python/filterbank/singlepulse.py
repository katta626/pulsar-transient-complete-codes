import csv

import numpy as np

from . import _lib as L
from .types import Header, Mask, Plan


def search(header, p, trials, start_sample=0, threshold=6.0, max_width=4096):
    t, ptr = L.f32(trials)
    params = {"threshold": threshold, "max_width": max_width}
    return L.take_json(L.lib.fb_sp_search(L.enc(Header(header).core()), L.enc(dict(p)), ptr, t.shape[0], t.shape[1], start_sample, L.enc(params)))


def run(paths, **cfg):
    if isinstance(paths, (str, bytes)) or hasattr(paths, "__fspath__"):
        paths = [paths]
    if "mask" in cfg and cfg["mask"] is not None:
        cfg["mask"] = dict(cfg["mask"])
    r = L.take_json(L.lib.fb_run_sp(L.enc([str(p) for p in paths]), L.enc(cfg)))
    return {"plan": Plan(r["plan"]), "mask": Mask(r["mask"]), "candidates": r["candidates"]}


def candidate_time(c, header):
    h = Header(header)
    return (c["epoch"]["imjd"] - h["start"]["imjd"]) * 86400.0 + (c["epoch"]["sec"] - h["start"]["sec"])


def to_csv(cands, path):
    if not cands:
        open(path, "w").close()
        return
    keys = ["dm", "snr", "sample", "width_samples", "width_s", "members", "dm_index", "beam"]
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["mjd", *keys])
        for c in cands:
            w.writerow([c["epoch"]["imjd"] + c["epoch"]["sec"] / 86400.0, *[c[k] for k in keys]])
