"""Independent spacetime gauge identity for the actual cloud/source bridge."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from environment_source import ThresholdCloud,angular_mode
from lorenz_ghp import KerrGHP
from lorenz_mode_jet import separated_jet
from lorenz_metric import spin1_metric,spin0_metric,nonstatic_metric


def test_actual_dipole_source_is_kg_of_lie_dragged_cloud_in_vacuum():
    cloud=ThresholdCloud(alpha=.3)
    for r in (10.,30.):
        theta=1.1
        g,one=spin1_metric(r,theta,20.,cloud.a,1,1,6,return_vector=True,full_current=True)
        _,zero=spin0_metric(r,theta,20.,cloud.a,1,1,6,return_vector=True)
        zeta=[-1j*(one[i]+zero[i])/g.omega for i in range(4)]
        gc=KerrGHP(r,theta,cloud.a,omega=cloud.omega,m=1,order=6)
        R,Rp=cloud.radial(r)[:,0]
        S,Sp,lam=angular_mode(theta,1,1,cloud.c2)
        eigenvalue=lam+cloud.a**2*cloud.omega**2-2*cloud.a*cloud.omega
        field=separated_jet(gc,0,eigenvalue,R,Rp,S,Sp,mass_squared=cloud.mu**2)
        dragged=sum(g.inv[i][j]*zeta[j]*gc.partial(field,i) for i in range(4) for j in range(4))
        total=KerrGHP(r,theta,cloud.a,omega=cloud.omega+g.omega,m=2,order=6)
        kg=(total.scalar_wave(dragged)-cloud.mu**2*dragged).value
        _,h=nonstatic_metric(r,theta,20.,cloud.a,1,1,6)
        source=cloud.lorenz_source(r,theta,[[v.value for v in row] for row in h])
        np.testing.assert_allclose(kg,source,rtol=2e-6,atol=1e-10)
