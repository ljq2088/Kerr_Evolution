#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from kerr_waveform_tools.background import isco_quantities


def main() -> None:
    parser = argparse.ArgumentParser(description="Write one validated lmax=6 BHPT ISCO-flux request")
    parser.add_argument("chi", type=float)
    parser.add_argument("output", type=Path)
    parser.add_argument("--flux-version", required=True)
    args = parser.parse_args()
    radius, _, _, _ = isco_quantities(args.chi)
    payload = {"chi": args.chi, "r_isco": radius, "lmax": 6, "flux_version": args.flux_version}
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="ascii")


if __name__ == "__main__":
    main()
