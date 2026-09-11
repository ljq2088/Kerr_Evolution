import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lorenz_static_boundary import StaticGreen


def test_static_green_recovers_regular_manufactured_solutions():
    for a,ell,kind in ((0.,0,'scalar'),(.6,2,'scalar'),(.877,4,'scalar'),(.6,1,'twoform'),(.6,2,'twoform')):
        def exact(r):
            t=r+.7;delta=r*r-2*r+a*a;dp=2*(r-1)
            if kind=='scalar':return np.array([1/t,-1/t**2,2/t**3])*(1+.3j)
            return np.array([delta/t**3,dp/t**3-3*delta/t**4,
                             2/t**3-6*dp/t**4+12*delta/t**5])*(1+.3j)
        def source(r):
            y,dy,ddy=exact(r)
            return (r*r-2*r+a*a)*ddy+(2*(r-1)*dy if kind=='scalar' else 0)-ell*(ell+1)*y
        green=StaticGreen(a,ell,source,kind=kind,pivot=6.,breaks=(4.,))
        for r in (green.rp+.02,3.5,15.,100.):
            value,_=green.evaluate(r)
            np.testing.assert_allclose(value,exact(r),rtol=2e-8,atol=3e-11,
                                       err_msg=f'{a=}, {ell=}, {kind=}, {r=}')


def test_regular_auxiliary_choice_preserves_vacuum_and_curvature():
    from lorenz_static_spin2 import sourced_static_spin2
    from lorenz_static_gauge import static_trace_metric
    from lorenz_tensor import trace,lorenz_constraint,linearized_einstein,extreme_weyl
    for a,ell,r in ((0.,2,4.5),(.6,2,9.),(.6,3,4.5)):
        g,h=sourced_static_spin2(r,1.1,a=a,ell=ell,boundary='regular')
        _,old=sourced_static_spin2(r,1.1,a=a,ell=ell)
        np.testing.assert_allclose(extreme_weyl(g,h),extreme_weyl(g,old),rtol=3e-10,atol=3e-12)
        np.testing.assert_allclose(trace(g,h).value,0,atol=3e-10)
        np.testing.assert_allclose([v.value for v in lorenz_constraint(g,h)],0,atol=3e-10)
        np.testing.assert_allclose([[v.value for v in row] for row in linearized_einstein(g,h)],0,atol=3e-10)
        if ell%2==0:
            g,h,tr=static_trace_metric(r,1.1,a=a,ell=ell,boundary='regular')
            np.testing.assert_allclose(trace(g,h).value,tr.value,rtol=3e-10,atol=3e-10)
            np.testing.assert_allclose([v.value for v in lorenz_constraint(g,h)],0,atol=3e-10)
            np.testing.assert_allclose([[v.value for v in row] for row in linearized_einstein(g,h)],0,atol=3e-10)


def test_free_scalar_has_prescribed_jumps_and_vacuum_metric():
    from lorenz_static_boundary import matched_free_scalar_metric
    from lorenz_static_gauge import _angular_jet
    from lorenz_tensor import trace,lorenz_constraint,linearized_einstein
    for a,ell in ((0.,2),(.6,4)):
        jumps=(.3+.1j,-.2+.05j)
        fields=[]
        for side in ('inside','outside'):
            g,h,field=matched_free_scalar_metric(6.,1.1,6.,a,ell,*jumps,side=side)
            fields.append(np.array([field.value,field.derivative(0).value]))
            np.testing.assert_allclose(trace(g,h).value,0,atol=3e-11)
            np.testing.assert_allclose([v.value for v in lorenz_constraint(g,h)],0,atol=3e-11)
            np.testing.assert_allclose([[v.value for v in row] for row in linearized_einstein(g,h)],0,atol=3e-11)
        np.testing.assert_allclose(fields[1]-fields[0],np.array(jumps)*_angular_jet(g,ell).value,atol=3e-14)


def test_regular_completion_quadrupole_preserves_vacuum_and_cancels_in_D_plus_E():
    from lorenz_completion import completion_metric
    from lorenz_tensor import trace,lorenz_constraint,linearized_einstein
    for r in (3.5,9.):
        differences=[]
        for mode in ('D','E'):
            g,h,tr=completion_metric(r,1.1,a=.6,mode=mode,y_boundary='regular-quadrupole')
            _,old,_=completion_metric(r,1.1,a=.6,mode=mode)
            differences.append(np.array([[h[i][j].value-old[i][j].value for j in range(4)] for i in range(4)]))
            np.testing.assert_allclose(trace(g,h).value,tr.value,atol=3e-10)
            np.testing.assert_allclose([v.value for v in lorenz_constraint(g,h)],0,atol=3e-10)
            np.testing.assert_allclose([[v.value for v in row] for row in linearized_einstein(g,h)],0,atol=3e-10)
        np.testing.assert_allclose(differences[0]+differences[1],0,atol=3e-13)


def test_dipole_completion_horizon_condition_from_fixed_horizon_geometry():
    from lorenz_completion import completion_metric
    for a in (0.,.6,.877):
        b=np.sqrt(1-a*a);rp=1+b
        maxima=[]
        for distance in (.01,.001):
            r=rp+distance;delta=r*r-2*r+a*a
            jac=np.eye(4);jac[0,1]=-(r*r+a*a)/delta;jac[3,1]=-a/delta
            h=np.zeros((4,4),complex)
            for mode,coefficient in (('F',1.),('G',2*b*b*rp),('C',-2*a*b*rp*rp)):
                _,piece,_=completion_metric(r,1.1,a=a,mode=mode,y_boundary='regular-quadrupole')
                h+=coefficient*np.array([[v.value for v in row] for row in piece])
            maxima.append(np.max(abs(jac.T@h@jac)))
        # A missed horizon pole grows by 10 or 100 under this refinement.
        assert maxima[1]<1.1*maxima[0]


def test_assembled_static_field_vacuum_and_completion_horizon():
    from environment_static_lorenz import StaticLorenzMode
    from lorenz_completion import completion_metric
    from lorenz_tensor import lorenz_constraint,linearized_einstein
    path=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/static_tetrad_a0.6_r6_L8_q20_j8_free8_paper.json'
    field=StaticLorenzMode(path)
    g,h=field.metric_jet(4.5,1.1)
    np.testing.assert_allclose([v.value for v in lorenz_constraint(g,h)],0,atol=2e-9)
    np.testing.assert_allclose([[v.value for v in row] for row in linearized_einstein(g,h)],0,atol=2e-9)
    maxima=[]
    for distance in (.01,.001):
        r=1.8+distance;delta=r*r-2*r+.36
        jac=np.eye(4);jac[0,1]=-(r*r+.36)/delta;jac[3,1]=-.6/delta
        total=np.zeros((4,4),complex)
        for mode,c in field.inside.items():
            _,h,_=completion_metric(r,1.1,a=.6,mode=mode,y_boundary='regular-quadrupole')
            total+=c*np.array([[v.value for v in row] for row in h])
        maxima.append(np.max(abs(jac.T@total@jac)))
    assert maxima[1]<1.1*maxima[0]
