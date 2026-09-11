import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from environment_angular_diagnostic import DenseRealHarmonic
from pybhpt.radial import RadialTeukolsky


def test_dense_real_harmonics_satisfy_angular_equation_and_cpp_separation():
    a=.8771530275949366;omega=2/(20**1.5+a);gamma=a*omega
    x,w=np.polynomial.legendre.leggauss(80);theta=np.arccos(x)
    for ell in (2,6,18):
        for spin in (-2,-1,0,1,2):
            h=DenseRealHarmonic(spin,ell,2,gamma)
            cpp=RadialTeukolsky(spin,ell,2,a,omega,np.array([30.]))
            assert abs(h.eigenvalue-cpp.eigenvalue)<5e-12
            s=h(theta);d=h(theta,deriv=1);dd=h(theta,deriv=2)
            A=h.eigenvalue-gamma**2+4*gamma
            potential=gamma**2*x*x-2*gamma*spin*x+spin+A-(2+spin*x)**2/(1-x*x)
            residual=dd+x/np.sqrt(1-x*x)*d+potential*s
            assert np.max(abs(residual))/max(1.,np.max(abs(potential*s)))<1e-11
            assert abs(2*np.pi*np.dot(w,abs(s)**2)-1)<1e-12
            repeated=DenseRealHarmonic(spin,ell,2,gamma)
            np.testing.assert_array_equal(h.couplingcoefficients,repeated.couplingcoefficients)


def test_dense_diagnostic_rejects_complex_frequency():
    with pytest.raises(ValueError,match='real frequency'):
        DenseRealHarmonic(2,2,2,.02+.01j)
