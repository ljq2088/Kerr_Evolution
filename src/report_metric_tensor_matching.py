"""All covariant metric components projected against smooth angular tests.

Quadratic local Taylor continuation estimates one-sided limits and a tenfold
smaller separation using the same independently computed vacuum jets.
Results remain diagnostics until angular and radial errors are converged.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import lpmv
from environment_source import kerr_metric
from lorenz_metric import nonstatic_metric


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--quadrature',type=int,default=12)
    parser.add_argument('--ellmax',type=int,default=8)
    parser.add_argument('--epsilon',type=float,default=5e-5)
    args=parser.parse_args()
    x,w=np.polynomial.legendre.leggauss(args.quadrature)
    weights=np.array([w*lpmv(2,j,x) for j in (2,3)])
    indices=[(a,b) for a in range(4) for b in range(a,4)]
    metric=kerr_metric(6.,np.pi/2,.6)
    op=1/(6**1.5+.6)
    ut=1/np.sqrt(-metric[0,0]-2*op*metric[0,3]-op*op*metric[3,3])
    ucov=metric@np.array([ut,0.,0.,op*ut])
    target=-8*(np.outer(ucov,ucov)+metric/2)/(ut*(36-12+.36))
    targets=np.array([[lpmv(2,j,0)*target[a,b] for a,b in indices] for j in (2,3)])
    # side, test function, radial Taylor derivative (0..2), tensor component
    sums=np.zeros((2,2,3,10),complex)
    cases=[]
    def encode(array):
        return np.stack((array.real,array.imag),axis=-1).tolist()
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'/f'metric_tensor_matching_q{args.quadrature}.json'
    for ell in range(2,args.ellmax+1):
        for side,r in enumerate((6.-args.epsilon,6.+args.epsilon)):
            for node,theta in enumerate(np.arccos(x)):
                _,h=nonstatic_metric(r,theta,6.,ell=ell,order=8)
                values=np.array([[h[a][b].value for a,b in indices],
                                 [h[a][b].derivative(0).value for a,b in indices],
                                 [h[a][b].derivative(0).derivative(0).value for a,b in indices]])
                sums[side]+=weights[:,node,None,None]*values[None,:,:]
        limits=[]
        for fraction in (1.,.1,0.):
            hv=[]
            dv=[]
            for side,sign in enumerate((-1.,1.)):
                d=sign*args.epsilon*(fraction-1.)
                hv.append(sums[side,:,0]+d*sums[side,:,1]+d*d*sums[side,:,2]/2)
                dv.append(sums[side,:,1]+d*sums[side,:,2])
            limits.append(dict(separation=2*args.epsilon*fraction,
                               value_jump=encode(hv[1]-hv[0]),
                               derivative_jump=encode(dv[1]-dv[0]),
                               derivative_error=encode(dv[1]-dv[0]-targets),
                               maximum_value_jump=float(np.max(abs(hv[1]-hv[0]))),
                               maximum_derivative_error=float(np.max(abs(dv[1]-dv[0]-targets)))))
        case=dict(ellmax=ell,limits=limits)
        cases.append(case)
        print(dict(ellmax=ell,maximum_limit_value_jump=limits[-1]['maximum_value_jump'],
                   maximum_limit_derivative_error=limits[-1]['maximum_derivative_error']),flush=True)
        out.write_text(json.dumps(dict(status='all_component_projected_diagnostic',jet_order=8,r0=6.,a=.6,m=2,
                 quadrature=args.quadrature,epsilon=args.epsilon,components=indices,tests=['P22','P32'],
                 expected_derivative_jump=targets.tolist(),cases=cases),indent=2)+'\n')


if __name__=='__main__':
    main()
