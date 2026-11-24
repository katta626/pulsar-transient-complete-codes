import json

import numpy as np

SECONDS_PER_DAY = 86400.0
K_DM = 1.0 / 2.41e-4


class Header(dict):
    @property
    def nchan(self):
        return len(self["freqs"]["mhz"])

    @property
    def npol(self):
        return self["npol"]

    @property
    def tsamp(self):
        return self["tsamp"]

    @property
    def nsamples(self):
        return self["nsamples"]

    @property
    def freqs(self):
        return np.asarray(self["freqs"]["mhz"])

    @property
    def foff(self):
        return self["freqs"]["chan_bw"]

    @property
    def mjd(self):
        return self["start"]["imjd"] + self["start"]["sec"] / SECONDS_PER_DAY

    @property
    def duration(self):
        return self["nsamples"] * self["tsamp"]

    @property
    def gaps(self):
        return self.get("gaps", [])

    def observation_id(self):
        return f"{self['telescope'].strip()}|{self['source_name'].strip()}|{self['start']['imjd']}|{self['start']['sec']:.3f}|beam{self['beam']}"

    def core(self):
        return {k: v for k, v in self.items() if k != "gaps"}

    def json(self):
        return json.dumps(self.core())

    def with_freqs(self, mhz, chan_bw=None):
        h = Header(self)
        mhz = [float(f) for f in mhz]
        bw = chan_bw if chan_bw is not None else (mhz[1] - mhz[0] if len(mhz) > 1 else self.foff)
        h["freqs"] = {"mhz": mhz, "chan_bw": bw}
        return h

    def intensity(self):
        h = Header(self)
        h["npol"] = 1
        h["pol_type"] = "Intensity"
        return h


class Mask(dict):
    @property
    def nchan(self):
        return len(self["chan_lo"])

    @property
    def nint(self):
        return self["nint"]

    @property
    def grid(self):
        return np.asarray(self["bad"], dtype=bool).reshape(self["nint"], self.nchan)

    @property
    def fraction(self):
        return float(np.mean(self["bad"])) if self["bad"] else 0.0

    def channel_fraction(self):
        return self.grid.mean(axis=0)

    def zap_freq_range(self, lo, hi):
        lo_e = np.asarray(self["chan_lo"])
        hi_e = np.asarray(self["chan_hi"])
        g = self.grid
        g[:, (hi_e > lo) & (lo_e < hi)] = True
        self["bad"] = g.reshape(-1).tolist()
        return self

    def union(self, other):
        if self["obs_id"] != other["obs_id"] or self.nchan != other.nchan or self.nint != other.nint:
            raise ValueError("masks do not share an observation and grid")
        m = Mask(self)
        m["bad"] = (self.grid | other.grid).reshape(-1).tolist()
        m["method"] = f"{self['method']}+{other['method']}"
        return m


class Plan(dict):
    @property
    def dms(self):
        return np.asarray(self["dms"])

    @property
    def ndm(self):
        return len(self["dms"])

    @property
    def f_ref(self):
        return self["f_ref"]
