import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lorenz_ghp import KerrGHP
from lorenz_corrector import corrector


def test_corrector_annihilates_conformal_test_tensor():
    # N is trace-free and anti-self-adjoint, hence N[g f]=0 for arbitrary f.
    g=KerrGHP(6.,1.2,.6,omega=.13,m=2,order=6)
    f=g.r*g.theta.sin()
    tensor=[[g.g[a][b]*f for b in range(4)] for a in range(4)]
    result=corrector(g,tensor)
    values=np.array([[x.value for x in row] for row in result])
    np.testing.assert_allclose(values,0,atol=2e-9)
