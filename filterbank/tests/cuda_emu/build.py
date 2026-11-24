import re
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def translate(src):
    src = re.sub(r"extern\s+__shared__\s+(\w+)\s+(\w+)\[\];", r"\1 *\2 = (\1 *)fbemu_dynamic_shared();", src)
    return re.sub(r"(\w+)<<<(.+?)>>>\(", r"fbemu_launch(fbemu_cfg(\2), \1, ", src)


def build(out_dir=None):
    out_dir = Path(out_dir or tempfile.mkdtemp())
    cpp = out_dir / "fbcuda_emu.cpp"
    cpp.write_text(translate((ROOT / "cuda" / "src" / "fbcuda.cu").read_text()))
    ext = "dylib" if sys.platform == "darwin" else "so"
    lib = out_dir / f"libfbcuda_emu.{ext}"
    cmd = ["c++", "-std=c++20", "-O1", "-shared", "-fPIC", "-I", str(HERE), "-I", str(ROOT / "cuda" / "include"), str(cpp), "-o", str(lib), "-lpthread"]
    subprocess.run(cmd, check=True)
    return lib


if __name__ == "__main__":
    print(build(sys.argv[1] if len(sys.argv) > 1 else None))
