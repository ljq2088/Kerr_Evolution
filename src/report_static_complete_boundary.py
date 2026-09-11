"""Boundary and local source diagnostics of the assembled static candidate."""
import json
from pathlib import Path
import numpy as np
from environment_static_lorenz import StaticLorenzMode
from report_static_tetrad_matching import tetrad_weights,spin_harmonic


def main():
    root=Path(__file__).resolve().parents[1]
    source=root/'docs/environment_reproduction/static_tetrad_a0.6_r6_L8_q20_j8_free8_paper.json'
    data=json.loads(source.read_text())
    mode=StaticLorenzMode(data)
    rp=1+np.sqrt(1-mode.a**2)
    boundary=[]
    for r in (rp+.01,rp+.001,30.,100.):
        h=mode(r,1.1);delta=r*r-2*r+mode.a**2
        jac=np.eye(4);jac[0,1]=-(r*r+mode.a**2)/delta;jac[3,1]=-mode.a/delta
        lengths=np.array([1,1,r,r])
        boundary.append(dict(r=r,advanced_maximum=float(np.max(abs(jac.T@h@jac))),
                             scaled_BL_maximum=float(np.max(abs(h/np.outer(lengths,lengths))))))
        print(boundary[-1],flush=True)
    x,w=np.polynomial.legendre.leggauss(12)
    pairs=[(i,j) for i in range(4) for j in range(i,4)]
    epsilon=5e-6
    jumps=[]
    for sign in (-1,1):
        r=mode.r0+sign*epsilon
        total=np.zeros((2,5,5),complex)
        for node,theta in enumerate(np.arccos(x)):
            _,h=mode.metric_jet(r,theta)
            values=np.array([h[i][j].value for i,j in pairs])
            derivatives=np.array([h[i][j].derivative(0).value for i,j in pairs])
            second=np.array([h[i][j].derivative(0).derivative(0).value for i,j in pairs])
            values+=-sign*epsilon*derivatives+epsilon**2*second/2
            derivatives-=sign*epsilon*second
            W,dW=tetrad_weights(mode.r0,theta,mode.a)
            tv=np.array([W@values,W@derivatives+dW@values])
            angular=np.array([[spin_harmonic(j,s,x[node])*w[node] for s in (0,2,1,0,0)] for j in range(5)])
            total+=tv[:,None,:]*angular[None,:,:]
        jumps.append(total)
        print('Completed assembled junction side',sign,flush=True)
    actual=jumps[1]-jumps[0]
    target=np.array(data['target'])[:,:5]
    target=target[...,0]+1j*target[...,1]
    r,a=mode.r0,mode.a
    scale=np.array([1,1/r**2,1/r**2,1/(r*r*(r*r-2*r+a*a)),1])[None,None,:]*np.array([1,r])[:,None,None]
    residual=abs((actual-target)*scale)
    result=dict(status='assembled_static_candidate_finite_resolution',matching_file=source.name,
                provenance=mode.provenance,boundary_samples=boundary,
                epsilon=epsilon,angular_nodes=12,test_degrees=list(range(5)),
                maximum_scaled_continuity_residual=float(np.max(residual[0])),
                maximum_scaled_derivative_residual=float(np.max(residual[1])),
                source_and_boundary_checks_do_not_establish_full_modal_convergence=True)
    out=root/'docs/environment_reproduction/static_complete_boundary.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k.startswith('maximum')}),flush=True)


if __name__=='__main__':
    main()
