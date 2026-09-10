from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .contracts import ContractError


def save_complex_series_plot(
    T,
    values,
    path: Path,
    *,
    symbol: str,
    title: str,
    event_lines: dict[str, float] | None = None,
    shaded_intervals: dict[str, tuple[float, float]] | None = None,
) -> Path:
    T = np.asarray(T, dtype=float)
    values = np.asarray(values, dtype=complex)
    if (T.ndim != 1 or values.shape != T.shape or len(T) < 2
            or not np.all(np.isfinite(T))
            or not np.all(np.isfinite(values.real))
            or not np.all(np.isfinite(values.imag))):
        raise ContractError("complex-series plotting requires finite one-dimensional arrays")
    figure, axis = plt.subplots(figsize=(10, 5.2), constrained_layout=True)
    axis.plot(T, values.real, color="#1769aa", lw=1.0, label=rf"$\mathrm{{Re}}\,{symbol}$")
    axis.plot(T, np.abs(values), color="#202020", lw=1.35, label=rf"$|{symbol}|$")
    for label, (start, stop) in (shaded_intervals or {}).items():
        axis.axvspan(start, stop, color="#6b7280", alpha=0.14, lw=0, label=label)
    for label, time in (event_lines or {}).items():
        if np.isfinite(time):
            axis.axvline(time, lw=0.8, ls=":", label=label)
    axis.axhline(0.0, color="#9ca3af", lw=0.6)
    axis.set(xlabel=r"$T/M$", ylabel="$" + symbol + "$", title=title)
    axis.set_xlim(0.0, float(T[-1]))
    axis.grid(alpha=0.18, lw=0.6)
    axis.legend(frameon=False, ncols=3)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


__all__ = ["save_complex_series_plot"]
