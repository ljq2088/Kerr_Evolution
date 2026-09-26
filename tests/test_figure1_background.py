import numpy as np
from types import SimpleNamespace
from li_field_green import LiFieldGreen
from environment_radial import RadialGreen
from report_figure1_background import field_inventory,field_grid

def test_inventory_retains_static_bound_and_negative_m():
    modes=field_inventory()
    assert len(modes)==88 and len(set(modes))==88
    assert len([(l,m) for l,m in modes if m==1])==5
    assert (12,-12) in modes and (12,12) in modes and (2,0) in modes

def test_grid_retains_orbit_boundary_and_integral_measure():
    c=SimpleNamespace(rp=1+np.sqrt(1-.88**2),a=.88,mu=.3,omega=.29629353472811226)
    for mg in (0,1,-1,13,-13):
        panels,r,w=field_grid(c,3.5,mg)
        assert 3.5 in panels and not np.any(r==3.5)
        assert np.all(np.diff(r)>0) and np.all(w>0)
        assert abs(sum(w)-(panels[-1]-panels[0]))<1e-10

def test_negative_frequency_green_matches_independent_retarded_solver():
    w=.29629353472811226-5/(3.5**1.5+.88)
    g=LiFieldGreen(.88,.3,w,4,-4,rmax=2000.)
    h=RadialGreen(.88,.3,w,4,-4,rmax=2000.,offset=1e-4,rtol=1e-11,infinity_method='coulomb')
    r=np.geomspace(2.,320.,100)
    K=g.insol.sol(3.)[0]*g.upsol.sol(r)[0]/g.w0
    H=h.insol.sol(3.)[0]*h.upsol.sol(r)[0]/h.w0
    assert np.linalg.norm(K-H)/np.linalg.norm(H)<1e-5
    assert np.max(abs(g.wronskian(r)/g.w0-1))<1e-6


def test_field_is_complex_coherent_sum_before_modulus():
    from render_figure1_background import evaluate
    from environment_source import angular_mode
    w=.45;theta=np.pi/2;c2=.88**2*(w*w-.3**2)
    s0=angular_mode(theta,2,0,c2)[0];s2=angular_mode(theta,2,2,c2)[0]
    modes={(2,0):(w,lambda r:np.ones_like(r,dtype=complex)*(-s2/s0)),
           (2,2):(w,lambda r:np.ones_like(r,dtype=complex))}
    assert abs(evaluate(modes,20.,theta,0.,2))<1e-13
    assert abs(evaluate(modes,20.,theta,np.pi/2,2))>1
