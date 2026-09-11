"""Explicit, process-local real-frequency angular backend for diagnostics.

This does not change the production backend on import. Install before any
source amplitudes are evaluated so source and vacuum factors use one solver.
"""
import numpy as np
from scipy.linalg import eigh
from pybhpt.swsh import SpinWeightedSpheroidalHarmonic


class DenseRealHarmonic(SpinWeightedSpheroidalHarmonic):
    def generate_eigs(self):
        if np.imag(self.spheroidicity)!=0:
            raise ValueError('This diagnostic requires real frequency')
        index=self.j-self.jmin
        size=index+round(20+abs(2*self.spheroidicity))+2
        matrix=self.sparse_matrix(size).toarray()
        if np.max(abs(matrix-matrix.T.conj()))>1e-12:
            raise ValueError('Angular matrix is not Hermitian')
        eigen,vectors=eigh(matrix,driver='evd')
        vector=vectors[:,index]*np.sign(vectors[index,index])
        return eigen[index],vector


def install_dense_angular_diagnostic():
    """Call once, before calculations, inside an isolated diagnostic process."""
    import lorenz_metric
    import lorenz_chi
    import lorenz_spin1
    import lorenz_spin1_chiral
    for module in (lorenz_metric,lorenz_chi,lorenz_spin1,lorenz_spin1_chiral):
        module.SpinWeightedSpheroidalHarmonic=DenseRealHarmonic
