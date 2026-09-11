"""Check auxiliary In/Up choices without claiming matched metric completion."""
import json
from pathlib import Path
import numpy as np
from lorenz_static_spin2 import sourced_static_spin2
from lorenz_static_gauge import static_trace_metric


def main():
    a,r0,theta,ell=.6,6.,1.1,2
    rp=1+np.sqrt(1-a*a)
    cases=[]
    for kind in ('spin2','trace'):
        for distance in (.01,.001,.0001):
            r=rp+distance
            if kind=='spin2':_,h=sourced_static_spin2(r,theta,r0=r0,a=a,ell=ell,boundary='regular')
            else:_,h,_=static_trace_metric(r,theta,r0=r0,a=a,ell=ell,boundary='regular')
            h=np.array([[v.value for v in row] for row in h])
            delta=r*r-2*r+a*a
            jac=np.eye(4);jac[0,1]=-(r*r+a*a)/delta;jac[3,1]=-a/delta
            advanced=jac.T@h@jac
            cases.append(dict(sector=kind,horizon_distance=distance,
                maximum_BL_component=float(np.max(abs(h))),maximum_advanced_component=float(np.max(abs(advanced))),
                advanced_components=advanced.real.tolist()))
    # The regular-minus-reference choice is a globally smooth homogeneous
    # gauge piece. Check both its value and radial derivative at the orbit.
    junction=[]
    epsilon=5e-5
    for kind in ('spin2','trace'):
        limits=[]
        for sign in (-1,1):
            r=r0+sign*epsilon
            pieces=[]
            for boundary in ('regular','reference'):
                if kind=='spin2':_,h=sourced_static_spin2(r,theta,r0=r0,a=a,ell=ell,boundary=boundary)
                else:_,h,_=static_trace_metric(r,theta,r0=r0,a=a,ell=ell,boundary=boundary,order=8)
                v=np.array([[v.value for v in row] for row in h])
                d=np.array([[v.derivative(0).value for v in row] for row in h])
                dd=np.array([[v.derivative(0).derivative(0).value for v in row] for row in h])
                shift=-sign*epsilon
                pieces.append(np.array([v+shift*d+shift*shift*dd/2,d+shift*dd]))
            limits.append(pieces[0]-pieces[1])
        jump=limits[1]-limits[0]
        junction.append(dict(sector=kind,maximum_homogeneous_difference_value_jump=float(np.max(abs(jump[0]))),
            maximum_homogeneous_difference_derivative_jump=float(np.max(abs(jump[1])))))
    result=dict(status='auxiliary_boundary_adaptation_only_completion_and_free_scalar_still_required',
        parameters=dict(a=a,r0=r0,theta=theta,ell=ell,epsilon=epsilon),
        horizon_cases=cases,junction_cases=junction,
        warning='Bounded sampled advanced components are evidence, not an analytic horizon limit or an error bound')
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/static_auxiliary_boundary.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(junction,indent=2),flush=True)


if __name__=='__main__':
    main()
