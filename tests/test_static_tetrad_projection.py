import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lorenz_ghp import KerrGHP
from report_static_tetrad_matching import tetrad_weights,spin_harmonic


def test_tetrad_projection_and_radial_product_rule_on_background():
    for a,r,theta in ((0.,4.,1.2),(.6,6.,.7),(.877,9.,1.5)):
        g=KerrGHP(r,theta,a,order=2)
        pairs=[(i,j) for i in range(4) for j in range(i,4)]
        h=np.array([g.g[i][j].value for i,j in pairs])
        dh=np.array([g.g[i][j].derivative(0).value for i,j in pairs])
        w,dw=tetrad_weights(r,theta,a)
        sigma=r*r+a*a*np.cos(theta)**2
        np.testing.assert_allclose(w@h,[0,0,0,2*sigma*sigma,4],atol=2e-10)
        np.testing.assert_allclose(dw@h+w@dh,[0,0,0,8*r*sigma,0],atol=2e-10)


def test_spin_harmonic_projection_normalization():
    x,w=np.polynomial.legendre.leggauss(32)
    for spin in (0,1,2):
        modes=np.array([spin_harmonic(ell,spin,x) for ell in range(spin,13)])
        np.testing.assert_allclose(2*np.pi*(modes*w)@modes.T,np.eye(len(modes)),atol=8e-14)
