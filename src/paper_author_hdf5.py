"""Read unmodified author HDF5 coefficients; no metric solver or guessed conventions.

Usage:
  mode = load_mode(2)
  radius, coeff = mode.nearest("In", 6.0)  # coeff[component, ell], ell=0..4
  radius, coeff = mode.at_index("Up", 0)  # exact outer-side orbit point

The first nine components use spin-weighted SPHERICAL harmonics. Writer code
exports q10 directly from the trace SPHEROIDAL radial expansion. No basis
conversion, interpolation, m substitution, or positive/negative-m summation is
performed here. ``contract`` requires caller-supplied angular values, including
an explicit q10 choice. Coefficients include ALL gravitational spheroidal modes
summed by the author calculation; the ell index on first nine is output spherical
ell, not an individual gravitational generating mode.
"""
from dataclasses import dataclass
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEPS = ROOT / "outputs/paper_metric_reference/python_deps"
if DEPS.exists():
    sys.path.insert(0, str(DEPS))
import h5py

DATA = ROOT / ("outputs/paper_metric_reference/ConorDyson_KerrLorenzMSF/"
    "GenerationCodes/Numerics-Asymptotics/h1-mmode-gen/h1mmodedat/"
    "data600/data600-80/data")
COMPONENTS = ("h_l+l+", "h_l-l-", "h_m+m+", "h_m-m-",
    "rho_h_l+m+", "rhob_h_l+m-", "rhob_h_l-m+", "rho_h_l-m-",
    "sigma_delta_h_l+l-", "h")
SPINS = (0, 0, 2, -2, 1, -1, 1, -1, 0, 0)

@dataclass(frozen=True)
class AuthorMode:
    path: Path
    m: int
    ell: np.ndarray
    radii: dict
    coefficients: dict

    def at_index(self, side, index):
        if side not in ("In", "Up"):
            raise ValueError("side must be In or Up")
        return float(self.radii[side][index]), self.coefficients[side][:, :, index].copy()

    def nearest(self, side, radius):
        rr = self.radii[side]
        if not min(rr) <= radius <= max(rr):
            raise ValueError(f"requested r={radius} outside {side} [{min(rr)}, {max(rr)}]")
        return self.at_index(side, int(np.argmin(abs(rr-radius))))

    def contract(self, side, index, angular_values):
        """Contract ten x ell coefficients with explicit ten x ell angular array.

        q10 must be supplied in its intended spheroidal basis; no default guess.
        Angular values include phi dependence if requested by caller. A single
        positive m complex Fourier coefficient is returned; no factor 2 is added.
        """
        radius, cc = self.at_index(side, index)
        angular_values = np.asarray(angular_values, dtype=complex)
        if angular_values.shape != cc.shape:
            raise ValueError(f"angular_values must have shape {cc.shape}")
        return radius, np.sum(cc * angular_values, axis=1)

def load_mode(m, path=None):
    path = Path(path) if path is not None else DATA / f"h1_a0.6_rp8.0_l4_m{m}.h5"
    rr, cc = {}, {}
    with h5py.File(path, "r") as f:
        group = f[f"m_{m}"]  # NEVER load another m for any component.
        for side in ("In", "Up"):
            rr[side] = np.array(group[f"r_{side.lower()}"], dtype=float)
            values = []
            for component in COMPONENTS:
                raw = np.asarray(group[side][component])
                if raw.dtype.names != ("Re", "Im"):
                    raise ValueError(f"unsupported complex storage: {raw.dtype}")
                values.append(raw["Re"] + 1j * raw["Im"])
            cc[side] = np.stack(values)
            if cc[side].ndim != 3 or cc[side].shape[2] != len(rr[side]):
                raise ValueError("unexpected component/radial dimensions")
            if not np.isfinite(cc[side]).all() or not np.isfinite(rr[side]).all():
                raise ValueError("nonfinite author data")
        if cc["In"].shape[:2] != cc["Up"].shape[:2]:
            raise ValueError("side harmonic sizes differ")
        if not (np.diff(rr["In"]) < 0).all() or not (np.diff(rr["Up"]) > 0).all():
            raise ValueError("unexpected side radial ordering")
    return AuthorMode(path, int(m), np.arange(cc["In"].shape[1]), rr, cc)
