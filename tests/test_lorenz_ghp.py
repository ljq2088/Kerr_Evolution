import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lorenz_jet import Jet
from lorenz_ghp import KerrGHP
from lorenz_mode_jet import separated_jet


def test_taylor_algebra_mixed_derivatives():
    r,t=Jet.variable(2.,0),Jet.variable(.7,1)
    f=(r**-1)*(2*t).sin()
    for i in range(4):
        for j in range(3):
            import math
            exact=(-1)**i*math.factorial(i)/2**(i+1)*2**j*np.sin(1.4+j*np.pi/2)
            np.testing.assert_allclose(f.derivative_value(i,j),exact,atol=1e-13)
    np.testing.assert_allclose(((r*r+1)/(r*r+1)).c[0,0],1)
    assert np.max(abs(((r*r+1)/(r*r+1)).c[1:]))<1e-12


def test_geometric_kinnersley_coefficients_and_scalar_operator():
    for radius,theta in ((3.5,.8),(20.,1.2)):
        ghp=KerrGHP(radius,theta,.6,omega=.2,m=2,order=5)
        tetra=np.array([[v.value for v in row] for row in ghp.tetrad])
        metric=np.array([[v.value for v in row] for row in ghp.g])
        eta=np.array([[0,-1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,1,0]])
        np.testing.assert_allclose(tetra@metric@tetra.T,eta,atol=1e-13)
        np.testing.assert_allclose(ghp.sc['eps'].value,0,atol=1e-13)
        np.testing.assert_allclose(ghp.sc['rho'].value,(-1/ghp.zeta).value,atol=1e-13)
        f=ghp.r**2*ghp.theta.sin()
        np.testing.assert_allclose(ghp.teukolsky(f,0).value,
                                   -.5*ghp.scalar_wave(f).value,rtol=1e-12,atol=1e-12)


def test_all_spin_operators_against_separated_equations():
    g=KerrGHP(6.,1.2,.6,omega=.13,m=2,order=6)
    for spin in (-2,-1,0,1,2):
        f=separated_jet(g,spin,6-spin*(spin+1),R0=1+.2j,R1=.3-.1j,S0=.7,S1=.2)
        residual=g.teukolsky(f,spin)
        # Check the equation and its first coordinate derivatives.
        for nr,nt in ((0,0),(1,0),(0,1)):
            np.testing.assert_allclose(residual.derivative_value(nr,nt),0,atol=2e-12,
                                       err_msg=f'spin {spin}, derivatives {nr},{nt}')
