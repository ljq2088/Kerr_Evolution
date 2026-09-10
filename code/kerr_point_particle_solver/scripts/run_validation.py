#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    subprocess.run([sys.executable, "-m", "pytest", "tests", "-q"],
                   cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
