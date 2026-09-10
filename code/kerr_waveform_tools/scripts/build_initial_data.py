#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from kerr_waveform_tools.transition import (
    FluxRecord,
    JsonInitialDataCache,
    TransitionSpec,
    generate_ori_thorne_initial_data,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build one fixed Ori--Thorne initial-data record from a versioned BHPT flux JSON")
    parser.add_argument("flux_json", type=Path)
    parser.add_argument("cache_directory", type=Path)
    args = parser.parse_args()
    record = FluxRecord(**json.loads(args.flux_json.read_text(encoding="ascii")))
    initial = generate_ori_thorne_initial_data(record)
    spec = TransitionSpec(record.chi, record.flux_version)
    key = spec.initial_data_key(initial.manifest.schema_version, initial.manifest.generator_version)
    JsonInitialDataCache(args.cache_directory, spec).store(key, initial)
    print(json.dumps({"key": key, "initial_data": initial.as_record()}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
