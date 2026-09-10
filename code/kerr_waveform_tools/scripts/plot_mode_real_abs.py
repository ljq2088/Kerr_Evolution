from __future__ import annotations

import argparse
from pathlib import Path

from kerr_waveform_tools.plotting import save_complex_series_plot
from kerr_waveform_tools.waveform_io import load_mode


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot the real part and magnitude of one saved psi4 mode.")
    parser.add_argument("mode_file", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--light-ring-scri", type=float)
    parser.add_argument("--source-off", type=float)
    args = parser.parse_args()

    T, psi4, metadata = load_mode(args.mode_file)
    ell, m = int(metadata["ell"]), int(metadata["m"])
    chi = metadata.get("background", {}).get("chi")
    title = rf"Peeling waveform: $(\ell,m)=({ell},{m})$"
    if chi is not None:
        title = rf"Peeling waveform: $\chi={chi:g}$, $(\ell,m)=({ell},{m})$"
    event_lines = {r"$T_{\rm stable}$": float(metadata["analysis_start"])}
    if args.light_ring_scri is not None:
        event_lines[r"$T_{\rm LR}^{\mathcal{I}^+}$"] = args.light_ring_scri
    if args.source_off is not None:
        event_lines[r"$T_{\rm source\ off}$"] = args.source_off
    save_complex_series_plot(
        T,
        psi4,
        args.output,
        symbol=rf"\psi_{{4,{ell}{m}}}",
        title=title,
        event_lines=event_lines,
    )


if __name__ == "__main__":
    main()
