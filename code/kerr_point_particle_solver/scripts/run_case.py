#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sminus2_point_particle.config import EvolutionConfig, OutputRequest
from sminus2_point_particle.initial_data import (
    FluxRecord,
    JsonInitialDataCache,
    generate_ori_thorne_initial_data,
)
from sminus2_point_particle.workflow import assemble_single_system, run_single_system
from sminus2_point_particle.run_storage import (
    data_run_paths,
    default_data_run_root,
    finalize_data_bundle,
    make_run_id,
)
from kerr_waveform_tools.transition import (
    TRANSITION_GENERATOR_VERSION,
    TRANSITION_SCHEMA_VERSION,
)


def repository_root() -> Path:
    return Path(__file__).resolve().parents[3]


def main() -> None:
    parser = argparse.ArgumentParser(description="Run one C1 fixed-m candidate system")
    parser.add_argument("config", type=Path)
    parser.add_argument("output_request", type=Path)
    parser.add_argument("flux", type=Path)
    parser.add_argument("run_directory", type=Path, nargs="?")
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--run-id")
    parser.add_argument("--initial-data-cache", type=Path)
    parser.add_argument("--t-end", type=float, help="Explicit short endpoint; omit to use T_LR^scri+post-LR")
    parser.add_argument("--restart", action="store_true")
    parser.add_argument("--no-progress", action="store_true", help="Disable the RK4 progress bar")
    args = parser.parse_args()
    if args.run_directory is not None and (args.output_root is not None or args.run_id is not None):
        parser.error("explicit run_directory cannot be combined with --output-root or --run-id")
    config = EvolutionConfig.from_json(args.config.read_text(encoding="ascii"))
    output = OutputRequest.from_json(args.output_request.read_text(encoding="ascii"))
    flux = FluxRecord(**json.loads(args.flux.read_text(encoding="ascii")))
    cache_directory = args.initial_data_cache or (
        repository_root() / "data" / "kerr_point_particle_evolution" / "initial_data"
    )
    key = config.initial_data_key(TRANSITION_SCHEMA_VERSION, TRANSITION_GENERATOR_VERSION)
    cache = JsonInitialDataCache(cache_directory, config)
    initial_data = cache.load(key)
    cache_status = "hit"
    if initial_data is None:
        cache_status = "miss-generated"
        initial_data = generate_ori_thorne_initial_data(flux)
        cache.store(key, initial_data)
    manifest = initial_data.manifest
    if (
        manifest.chi != flux.chi
        or manifest.lmax != flux.lmax
        or manifest.flux_model_version != flux.flux_version
        or manifest.flux_total_over_mu2 != flux.flux_total_over_mu2
        or manifest.relative_lmax_increment != flux.relative_lmax_increment
        or manifest.flux_provenance != flux.provenance
    ):
        raise ValueError("initial-data cache record does not match the supplied flux record")
    cache_path = cache_directory / f"{key}.json"
    provenance = {
        "key": key,
        "status": cache_status,
        "record_path": str(cache_path.resolve()),
        "record_sha256": hashlib.sha256(cache_path.read_bytes()).hexdigest(),
        "schema_version": initial_data.manifest.schema_version,
        "generator_version": initial_data.manifest.generator_version,
        "flux_model_version": initial_data.manifest.flux_model_version,
    }
    system = assemble_single_system(
        config,
        initial_data,
        output,
        initial_data_provenance=provenance,
    )
    managed_paths = None
    if args.run_directory is None:
        run_id = args.run_id or make_run_id(config)
        managed_paths = data_run_paths(args.output_root or default_data_run_root(), run_id)
        run_directory = managed_paths.staging
    else:
        run_directory = args.run_directory
    result = run_single_system(system, run_directory, T_end=args.t_end, restart=args.restart,
                               show_progress=not args.no_progress)
    if managed_paths is not None and args.t_end is None:
        final_directory = finalize_data_bundle(managed_paths, config, output)
        result = type(result)(final_directory, result.final_time, result.final_state, result.output_count)
    print(json.dumps({"directory": str(result.directory), "final_time": result.final_time,
                      "output_count": result.output_count}, sort_keys=True))


if __name__ == "__main__":
    main()
