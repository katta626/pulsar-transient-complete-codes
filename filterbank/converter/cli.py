from __future__ import annotations

import argparse
from typing import Optional, Sequence

from .converter import ConversionError, convert_psrfits_to_filterbank


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(prog="fits2fil", description="Convert PSRFITS search data (or filterbank) to SIGPROC filterbank.")
    p.add_argument("inputs", nargs="+")
    p.add_argument("-o", "--output")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--nbits", type=int)
    p.add_argument("--scales", choices=["auto", "always", "never"], default="auto")
    p.add_argument("--zero-off", choices=["spec", "dspsr", "ignore"])
    p.add_argument("--keep-pols", action="store_true")
    p.add_argument("--order", choices=["descending", "keep"], default="descending")
    a = p.parse_args(argv)
    try:
        out = convert_psrfits_to_filterbank(
            a.inputs, a.output, overwrite=a.overwrite, nbits=a.nbits, scales=a.scales,
            zero_off=a.zero_off, keep_pols=a.keep_pols, order=a.order,
        )
    except ConversionError as exc:
        p.exit(status=2, message=f"fits2fil: error: {exc}\n")
    print(f"Wrote {out}")
    return 0
