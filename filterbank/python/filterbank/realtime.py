import time
from pathlib import Path

import numpy as np

from . import dedisp, rfi, singlepulse
from .io import Reader


def follow(path, gulp=65536, overlap=0, poll=0.5, idle_timeout=30.0, **open_kw):
    path = Path(path)
    pos = 0
    last_growth = time.time()
    while True:
        r = Reader(path, **open_kw)
        n = r.header.nsamples
        grew = False
        while n >= pos + gulp + overlap:
            yield pos, r.read(pos, gulp + overlap)
            pos += gulp
            grew = True
        r.close()
        if grew:
            last_growth = time.time()
        elif time.time() - last_growth > idle_timeout:
            tail = None
            r = Reader(path, **open_kw)
            if r.header.nsamples > pos + overlap:
                tail = r.read(pos, r.header.nsamples - pos)
            r.close()
            if tail is not None and len(tail):
                yield pos, tail
            return
        time.sleep(poll)


def watch_directory(directory, pattern=("*.fil", "*.dada"), poll=2.0, idle_timeout=None):
    seen = set()
    start = time.time()
    while True:
        new = []
        for pat in pattern:
            for p in sorted(Path(directory).glob(pat)):
                if p not in seen:
                    seen.add(p)
                    new.append(p)
        for p in new:
            yield p
        if idle_timeout is not None and not new and time.time() - start > idle_timeout:
            return
        time.sleep(poll)


class Searcher:
    def __init__(self, header, dm_lo=0.0, dm_hi=1000.0, dm_step=None, threshold=7.0, max_width=1024, backend="auto", interval=1.0, on_candidate=None):
        self.h = header
        self.plan = dedisp.plan(header, dm_lo, dm_hi, step=dm_step)
        self.table, self.maxd = dedisp.delay_table(self.plan)
        self.threshold = threshold
        self.max_width = max_width
        self.backend = backend
        self.interval = interval
        self.on_candidate = on_candidate or (lambda c: None)
        self.mask = rfi.new_mask(header, interval)
        self.ncands = 0

    def process(self, start, block):
        data = block[:, 0, :] if block.ndim == 3 else block
        if data.shape[0] <= self.maxd + 16:
            return []
        self.mask = rfi.find(self.h, data, start_sample=start, mask=self.mask, interval=self.interval)
        cleaned, _ = rfi.apply(self.h, data, self.mask, start_sample=start)
        kill = self.mask.channel_fraction() > 0.5
        trials = dedisp.dedisperse(self.h, self.plan, cleaned, killmask=kill, backend=self.backend)
        cands = singlepulse.search(self.h, self.plan, trials, start_sample=start, threshold=self.threshold, max_width=self.max_width)
        for c in cands:
            self.on_candidate(c)
        self.ncands += len(cands)
        return cands


def run_live(path, gulp=65536, poll=0.5, idle_timeout=30.0, **kw):
    with Reader(path) as r:
        header = r.header
    searcher = Searcher(header, **kw)
    allc = []
    for start, block in follow(path, gulp=gulp, overlap=searcher.maxd, poll=poll, idle_timeout=idle_timeout):
        allc.extend(searcher.process(start, block))
    return allc
