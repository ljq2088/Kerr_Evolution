#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np

from sminus2_point_particle.strain import fixed_frequency_integrate, write_ffi_result
from sminus2_point_particle.io import load_mode, load_run_metadata


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert uniform peeling psi4_lm samples to conditional H_lm=r*h_lm using FFI")
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--ell", type=int, required=True)
    parser.add_argument("--m", type=int, required=True)
    parser.add_argument("--phi-dot-at-start", type=float, help="Optional cross-check against mode metadata")
    parser.add_argument("--test-input", action="store_true", help="Allow a mode explicitly marked test; output remains test")
    parser.add_argument("--cutoff-factor", type=float, default=0.75)
    parser.add_argument("--taper-width", type=float, default=10.0)
    args = parser.parse_args()
    T, psi4, source_metadata = load_mode(args.input, expected_ell=args.ell, expected_m=args.m,
                                         allow_incomplete=args.test_input)
    source_metadata = dict(source_metadata)
    if "background" not in source_metadata or "amplitude_scaling" not in source_metadata:
        run_metadata = load_run_metadata(args.input.parent, allow_incomplete=args.test_input)
        if run_metadata.get("config_hash") != source_metadata["config_hash"]:
            raise ValueError("sibling run metadata does not match the mode config hash")
        config = run_metadata["config"]
        source_metadata["background"] = {
            "mass_scale": config["mass_scale"],
            "M_internal": 1.0,
            "chi": config["chi"],
            "a_internal": config["chi"],
            "L_internal": 1.0,
        }
        source_metadata["amplitude_scaling"] = (
            "unit particle-mass response; no physical mu, luminosity-distance, "
            "or BBH amplitude mapping applied"
        )
    if source_metadata["status"] != "complete" and not (args.test_input and source_metadata["status"] == "test"):
        raise ValueError("converter requires a complete full run unless --test-input is explicit")
    if source_metadata["run_kind"] != ("test" if args.test_input else "full"):
        raise ValueError("mode run_kind is inconsistent with converter mode")
    phi_dot = float(source_metadata["phi_dot_at_analysis_start"])
    analysis_start = float(source_metadata["analysis_start"])
    if not T[0] <= analysis_start < T[-1]:
        raise ValueError("mode data do not cover the recorded C2 analysis start")
    if args.phi_dot_at_start is not None and not np.isclose(args.phi_dot_at_start, phi_dot, rtol=1e-12, atol=1e-14):
        raise ValueError("CLI phi-dot cross-check does not match mode metadata")
    if not np.isclose(T[-1], source_metadata["T1"], rtol=0.0, atol=1e-12):
        raise ValueError("mode data do not terminate at recorded T1")
    result = fixed_frequency_integrate(T, psi4, ell=args.ell, m=args.m,
                                       phi_dot_at_start=phi_dot, transform_start=float(T[0]),
                                       analysis_start=analysis_start,
                                       stop=float(source_metadata["T1"]), input_schema=source_metadata["schema"],
                                       cutoff_factor=args.cutoff_factor, taper_width=args.taper_width,
                                       source_input_hash=source_metadata["input_hash"],
                                       source_events=source_metadata["events"],
                                       source_normalization=source_metadata["normalization"],
                                       source_error_metadata=source_metadata["error_metadata"],
                                       source_environment=source_metadata["environment"],
                                       source_background=source_metadata["background"],
                                       amplitude_scaling=source_metadata["amplitude_scaling"],
                                       sampling_status="owner-waived")
    write_ffi_result(args.output, result, status="test" if args.test_input else "complete",
                     source_mode_metadata=source_metadata)


if __name__ == "__main__":
    main()
