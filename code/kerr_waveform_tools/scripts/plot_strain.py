#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from kerr_waveform_tools.strain import load_ffi_result


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot the real part and magnitude of one C2 H_lm result.")
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--analysis-only", action="store_true")
    args = parser.parse_args()

    result, metadata = load_ffi_result(args.input)
    T, H = result.T, result.H
    if args.analysis_only:
        selected = (T >= metadata["analysis_start"]) & (T <= metadata["analysis_stop"])
        T, H = T[selected], H[selected]

    figure, axis = plt.subplots(figsize=(10, 5.2), constrained_layout=True)
    axis.plot(T, H.real, color="#1769aa", lw=1.05, label=rf"$\mathrm{{Re}}\,H_{{{metadata['ell']}{metadata['m']}}}$")
    axis.plot(T, np.abs(H), color="#202020", lw=1.4, label=rf"$|H_{{{metadata['ell']}{metadata['m']}}}|$")
    if not args.analysis_only:
        axis.axvspan(0.0, metadata["analysis_start"], color="#9ca3af", alpha=0.12, lw=0,
                     label="outside analysis interval")
    axis.axvspan(metadata["transform_start"], metadata["transform_start"] + metadata["taper_width"],
                 color="#6b7280", alpha=0.20, lw=0)
    axis.axvspan(metadata["recommended_stop"], metadata["analysis_stop"],
                 color="#6b7280", alpha=0.16, lw=0, label="tapered edge")
    axis.axhline(0.0, color="#9ca3af", lw=0.6)
    chi = metadata["source_background"]["chi"]
    axis.set(
        xlabel=r"$T/M$",
        ylabel=rf"$H_{{{metadata['ell']}{metadata['m']}}}=r h_{{{metadata['ell']}{metadata['m']}}}$",
        title=rf"C2 FFI strain: $\chi={chi:g}$, $(\ell,m)=({metadata['ell']},{metadata['m']})$",
    )
    axis.legend(frameon=False, ncols=3)
    axis.grid(alpha=0.18, lw=0.6)
    if not args.analysis_only:
        axis.set_xlim(0.0, metadata["analysis_stop"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=180)
    plt.close(figure)


if __name__ == "__main__":
    main()
