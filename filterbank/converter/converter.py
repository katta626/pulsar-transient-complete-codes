from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional, Sequence, Union

ROOT = Path(__file__).resolve().parents[1]


class ConversionError(RuntimeError):
    pass


def find_fb() -> str:
    env = os.environ.get("FB_BIN")
    if env:
        return env
    for p in (ROOT / "target" / "release" / "fb", ROOT / "target" / "debug" / "fb"):
        if p.exists():
            return str(p)
    found = shutil.which("fb")
    if found:
        return found
    raise ConversionError("fb binary not found; run `cargo build --release` in the filterbank directory")


def convert_psrfits_to_filterbank(
    input_paths: Sequence[Union[str, Path]],
    output_path: Optional[Union[str, Path]] = None,
    *,
    overwrite: bool = False,
    nbits: Optional[int] = None,
    scales: str = "auto",
    zero_off: Optional[str] = None,
    keep_pols: bool = False,
    order: str = "descending",
) -> Path:
    paths = [Path(p).expanduser().resolve() for p in input_paths]
    out = Path(output_path).expanduser().resolve() if output_path else paths[0].with_suffix(".fil")
    cmd = [find_fb(), "convert", *map(str, paths), "-o", str(out), "--scales", scales, "--order", order]
    if overwrite:
        cmd.append("--overwrite")
    if nbits is not None:
        cmd += ["--nbits", str(nbits)]
    if zero_off is not None:
        cmd += ["--zero-off", zero_off]
    if keep_pols:
        cmd.append("--keep-pols")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise ConversionError(r.stderr.strip())
    return out
