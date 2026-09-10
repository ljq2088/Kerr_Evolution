#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from sminus2_point_particle.config import EvolutionConfig, OutputRequest
from sminus2_point_particle.evolution import EvolutionState, SpatialRHS
from sminus2_point_particle.initial_data import FluxRecord, generate_ori_thorne_initial_data
from sminus2_point_particle.workflow import assemble_single_system


def median_call_time(function, repeats: int, warmup: int) -> float:
    for _ in range(warmup):
        function()
    samples = []
    for _ in range(5):
        start = perf_counter()
        for _ in range(repeats):
            function()
        samples.append((perf_counter() - start) / repeats)
    return float(np.median(samples))


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark the fixed-m RHS hot path.")
    parser.add_argument("config", type=Path)
    parser.add_argument("output_request", type=Path)
    parser.add_argument("flux", type=Path)
    parser.add_argument("--repeats", type=int, default=50)
    parser.add_argument("--warmup", type=int, default=3)
    parser.add_argument("--n-r", type=int)
    parser.add_argument("--n-y", type=int)
    args = parser.parse_args()
    if args.repeats <= 0 or args.warmup < 0:
        parser.error("repeats must be positive and warmup nonnegative")

    config = EvolutionConfig.from_json(args.config.read_text(encoding="ascii"))
    if args.n_r is not None or args.n_y is not None:
        config = replace(config, n_r=args.n_r or config.n_r, n_y=args.n_y or config.n_y)
    output = OutputRequest.from_json(args.output_request.read_text(encoding="ascii"))
    flux = FluxRecord(**json.loads(args.flux.read_text(encoding="ascii")))
    system = assemble_single_system(config, generate_ori_thorne_initial_data(flux), output)

    profile_R = np.sin(3 * system.grid.RR) * (1 - system.grid.yy**2)
    profile_P = np.cos(2 * system.grid.RR) * (1 - system.grid.yy**2)
    state = EvolutionState(profile_P.astype(complex), profile_R.astype(complex))
    field_rhs = SpatialRHS(system.grid, config.chi, config.m, source=None)
    T_on = max(config.tau_on_over_m, system.timestep.t_stable)
    T_off = np.nextafter(system.trajectory.events.source_off, np.inf)

    timings = {
        "trajectory_sample": median_call_time(lambda: system.trajectory.sample(T_on), args.repeats, args.warmup),
        "source_total_on": median_call_time(lambda: system.source(T_on), args.repeats, args.warmup),
        "field_rhs": median_call_time(lambda: field_rhs(T_on, state), args.repeats, args.warmup),
        "full_rhs_on": median_call_time(lambda: system.rhs(T_on, state), args.repeats, args.warmup),
        "full_rhs_off": median_call_time(lambda: system.rhs(T_off, state), args.repeats, args.warmup),
    }
    record = {
        "chi": config.chi,
        "m": config.m,
        "shape": list(system.grid.RR.shape),
        "T_on": T_on,
        "T_off": T_off,
        "repeats_per_sample": args.repeats,
        "uses_local_support": system.source.uses_local_support,
        "median_seconds_per_call": timings,
        "source_fraction_of_on_rhs": timings["source_total_on"] / timings["full_rhs_on"],
    }
    print(json.dumps(record, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
