import ctypes
import os
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
_EXT = {"darwin": "dylib", "win32": "dll"}.get(sys.platform, "so")


def cuda_library():
    for c in (os.environ.get("FB_CUDA_LIB"), _ROOT / "cuda" / "build" / f"libfbcuda.{_EXT}", Path(__file__).parent / f"libfbcuda.{_EXT}"):
        if c and Path(c).exists():
            return str(c)
    return None


_cuda = None


def cuda():
    global _cuda
    if _cuda is None:
        path = cuda_library()
        if path is None:
            return None
        try:
            lib = ctypes.CDLL(path)
            lib.fbcuda_device_count.restype = ctypes.c_int
            if lib.fbcuda_device_count() <= 0:
                return None
            _cuda = lib
        except OSError:
            return None
    return _cuda


def available():
    out = ["cpu"]
    if cuda() is not None:
        out.append("cuda")
    return out


def resolve(backend="auto"):
    if backend == "auto":
        return "cuda" if cuda() is not None else "cpu"
    if backend not in available():
        raise RuntimeError(f"backend '{backend}' is not available here (have: {', '.join(available())})")
    return backend
