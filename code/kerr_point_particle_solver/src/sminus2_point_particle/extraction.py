from __future__ import annotations

from typing import Callable
import numpy as np

from kerr_waveform_tools.harmonics import (
    project_mode,
    projection_row,
    spin_weighted_spherical_harmonic,
)

from .errors import ContractError, OpenBoundaryError


def scri_slice(
    field: np.ndarray,
    radial_nodes=None,
    radial_axis: int = 0,
    endpoint_reconstructor: Callable[[np.ndarray, np.ndarray, int], np.ndarray] | None = None,
) -> np.ndarray:
    field = np.asarray(field)
    if radial_nodes is None:
        raise ContractError("radial_nodes are required to identify R=0 without ambiguity")
    radial_nodes = np.asarray(radial_nodes, dtype=float)
    if field.shape[radial_axis] != radial_nodes.size:
        raise ContractError("radial field axis does not match radial_nodes")
    matches = np.flatnonzero(np.isclose(radial_nodes, 0.0, rtol=0.0, atol=1e-15))
    if matches.size:
        return np.take(field, int(matches[0]), axis=radial_axis)
    if endpoint_reconstructor is None:
        raise OpenBoundaryError("this grid needs an explicit SCRI+ endpoint reconstructor")
    return endpoint_reconstructor(radial_nodes, field, radial_axis)


def simpson_projection_row(grid, ell: int, m: int) -> np.ndarray:
    from .grid import simpson_weights_augmented
    nodes, weights = simpson_weights_augmented(grid)
    return projection_row(ell, m, nodes, weights)


def project_scri_mode(field: np.ndarray, grid, ell: int, m: int) -> complex:
    from .grid import reconstruct_axis_endpoints
    scri = scri_slice(field, grid.R)
    augmented = reconstruct_axis_endpoints(scri, grid)
    return project_mode(augmented, simpson_projection_row(grid, ell, m))
