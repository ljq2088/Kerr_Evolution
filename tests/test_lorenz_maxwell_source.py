import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lorenz_ghp import KerrGHP
from lorenz_tensor import vector_covariant_derivative,tensor_divergence


def test_maxwell_decoupling_source_normalization():
    g=KerrGHP(8.,1.1,.6,omega=.13,m=2,order=6)
    xi=[g.r*g.theta.sin(),g.r**2,g.theta.cos(),g.r*g.theta.cos()]
    d=vector_covariant_derivative(g,xi)
    F=[[d[b][a]-d[a][b] for b in range(4)] for a in range(4)]
    j=tensor_divergence(g,F)
    l,n,m,b=g.tetrad
    def proj(u):
        return sum(u[a]*j[a] for a in range(4))
    jl,jn,jm,jb=map(proj,(l,n,m,b))
    sc=g.sc
    source0=(g.derivative(jl,1,1,'eth')[0]-(2*sc['tau']+sc['taup'].conjugate())*jl
             -g.derivative(jm,1,-1,'thorn')[0]+(2*sc['rho']+sc['rho'].conjugate())*jm)/2
    source2=(-g.derivative(jn,-1,-1,'ethp')[0]+(2*sc['taup']+sc['tau'].conjugate())*jn
             +g.derivative(jb,-1,1,'thornp')[0]-(2*sc['rhop']+sc['rhop'].conjugate())*jb)/2
    for s,u,v,source in ((1,l,m,source0),(-1,b,n,source2)):
        phi=sum(F[a][c]*u[a]*v[c] for a in range(4) for c in range(4))
        np.testing.assert_allclose(g.teukolsky(phi,s).value,source.value,rtol=1e-10,atol=1e-11)
