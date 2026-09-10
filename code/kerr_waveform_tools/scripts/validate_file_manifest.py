#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate the files table in a research-data manifest.")
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="ascii"))
    root = args.manifest.parent
    failures = []
    for relative, expected in manifest["files"].items():
        path = root / relative
        actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
        if actual != expected:
            failures.append({"path": relative, "expected": expected, "actual": actual})
    if failures:
        raise SystemExit(json.dumps({"status": "failed", "failures": failures}, indent=2))
    print(json.dumps({"status": "pass", "files": len(manifest["files"])}, sort_keys=True))


if __name__ == "__main__":
    main()
