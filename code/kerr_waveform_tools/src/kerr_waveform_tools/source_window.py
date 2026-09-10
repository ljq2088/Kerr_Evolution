from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .contracts import ContractError


@dataclass(frozen=True)
class SourceWindowJet:
    value: float
    first: float
    second: float


def quintic_source_turn_on(T: float, tau_on: float = 20.0) -> SourceWindowJet:
    T = float(T)
    tau_on = float(tau_on)
    if not np.isfinite(T) or not np.isfinite(tau_on) or not 0 < tau_on <= 40:
        raise ContractError("source turn-on requires finite T and 0 < tau_on <= 40M")
    if T <= 0:
        return SourceWindowJet(0.0, 0.0, 0.0)
    if T >= tau_on:
        return SourceWindowJet(1.0, 0.0, 0.0)
    x = T / tau_on
    return SourceWindowJet(
        10 * x**3 - 15 * x**4 + 6 * x**5,
        30 * x**2 * (1 - x) ** 2 / tau_on,
        60 * x * (1 - x) * (1 - 2 * x) / tau_on**2,
    )


def quintic_source_turn_on_values(T, tau_on: float = 20.0) -> np.ndarray:
    values = np.asarray(T, dtype=float)
    if (
        not np.all(np.isfinite(values))
        or not np.isfinite(tau_on)
        or not 0 < tau_on <= 40
    ):
        raise ContractError("source turn-on requires finite T and 0 < tau_on <= 40M")
    x = np.clip(values / float(tau_on), 0.0, 1.0)
    return np.asarray(10 * x**3 - 15 * x**4 + 6 * x**5, dtype=float)


__all__ = [
    "SourceWindowJet",
    "quintic_source_turn_on",
    "quintic_source_turn_on_values",
]
