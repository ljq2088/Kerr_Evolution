"""Build an isolated patched nu solver, linked to the unchanged installed MST kernel.

No pip install, LD_PRELOAD, or replacement of the production extension occurs.
Run this script and the validation report in separate processes. The shared
library uses -Bsymbolic for its patched nu implementation and is loaded LOCAL.
"""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import urllib.request
import cybhpt_full

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "outputs/pybhpt_mst_validation"
COMMIT = "9fe9c57e2d1c92d944ba70ba3c1b81b665e9d274"
SOURCE_SHA = "0688632c4d54e4f4faa1fc5267c7a35bce51ae524b86b6c02bfa35d5d53dbb28"

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def build():
    BASE.mkdir(parents=True, exist_ok=True)
    inc = BASE / "include"
    inc.mkdir(exist_ok=True)
    pending = ["mst.hpp", "nusolver.hpp"]
    seen = set()
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        seen.add(name)
        dest = inc / name
        if not dest.exists():
            url = f"https://raw.githubusercontent.com/znasipak/pybhpt/{COMMIT}/cpp/include/{name}"
            with urllib.request.urlopen(url, timeout=30) as response:
                dest.write_bytes(response.read())
        pending.extend(x for x in re.findall(r'#include\s+"([^"]+)"', dest.read_text()) if "/" not in x)
    original = ROOT / "outputs/lorenz_reference/pybhpt_radial_source/v1.0.0/nusolver.cpp"
    if not original.exists():
        original.parent.mkdir(parents=True, exist_ok=True)
        url = f"https://raw.githubusercontent.com/znasipak/pybhpt/{COMMIT}/cpp/src/nusolver.cpp"
        with urllib.request.urlopen(url, timeout=30) as response:
            original.write_bytes(response.read())
    if sha(original) != SOURCE_SHA:
        raise RuntimeError("Reference nusolver.cpp differs from the audited v1.0.0 source")
    patched = BASE / "nusolver_patched.cpp"
    patched.write_bytes(original.read_bytes())
    patch = ROOT / "docs/environment_reproduction/patches/pybhpt_v1_nu_low_frequency.patch"
    subprocess.run(["patch", "--batch", str(patched), str(patch)], check=True)
    bridge = ROOT / "src/pybhpt_mst_validation_bridge.cpp"
    extension = Path(cybhpt_full.__file__).resolve()
    extension_before = sha(extension)
    output = BASE / "libpatched_mst.so"
    command = ["g++", "-shared", "-fPIC", "-O2", "-std=c++14", "-Wl,-Bsymbolic",
               "-I" + str(inc), str(patched), str(bridge), str(extension), "-o", str(output)]
    subprocess.run(command, check=True)
    assert sha(extension) == extension_before
    manifest = {"reference_commit": COMMIT, "builder_script_sha256": sha(__file__), "original_nusolver_sha256": SOURCE_SHA,
                "patched_nusolver_sha256": sha(patched), "patch_sha256": sha(patch),
                "bridge_sha256": sha(bridge), "library_sha256": sha(output),
                "installed_extension_path": str(extension), "installed_extension_sha256": extension_before,
                "headers": {name: sha(inc/name) for name in sorted(seen)}, "command": command,
                "scope": "isolated patched nu translation unit plus unchanged installed MST series kernel"}
    (BASE / "build_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(output)

if __name__ == "__main__":
    build()
