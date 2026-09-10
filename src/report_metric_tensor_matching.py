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
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--smooth',action='store_true',help='Use sin(theta) d_theta, a globally smooth angular vector')
    args=parser.parse_args()
    x,w=np.polynomial.legendre.leggauss(args.quadrature)
    weights=np.array([w*lpmv(2,j,x) for j in (2,3)])
    indices=[(a,b) for a in range(4) for b in range(a,4)]
    theta_count=np.array([int(a==2)+int(b==2) for a,b in indices])
    polar_factor=(np.sqrt(1-x*x)[:,None]**theta_count[None,:]
                  if args.smooth else np.ones((len(x),10)))
    weights=weights[:,:,None]*polar_factor[None,:,:]
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
    suffix='_smooth' if args.smooth else ''
    if args.epsilon!=5e-5:
        suffix+=f'_eps{args.epsilon:g}'
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'/f'metric_tensor_matching_q{args.quadrature}{suffix}.json'
    baseline=None
    first_ell=2
    if args.resume:
        previous=json.loads(out.read_text())
        expected=dict(jet_order=8,r0=6.,a=.6,m=2,quadrature=args.quadrature,epsilon=args.epsilon)
        if any(previous.get(k)!=v for k,v in expected.items()):
            raise ValueError('Resume parameters do not match the saved calculation')
        if previous.get('smooth_polar_tests',False)!=args.smooth:
            raise ValueError('Resume angular test vectors differ')
        np.testing.assert_allclose(previous['expected_derivative_jump'],targets,rtol=1e-14,atol=1e-14)
        cases=previous['cases']
        if not cases:
            raise ValueError('No completed mode to resume')
        baseline=cases[-1]['limits']
        first_ell=cases[-1]['ellmax']+1
    def decode(array):
        values=np.asarray(array)
        return values[...,0]+1j*values[...,1]
    for ell in range(first_ell,args.ellmax+1):
        for side,r in enumerate((6.-args.epsilon,6.+args.epsilon)):
            for node,theta in enumerate(np.arccos(x)):
                _,h=nonstatic_metric(r,theta,6.,ell=ell,order=8)
                values=np.array([[h[a][b].value for a,b in indices],
                                 [h[a][b].derivative(0).value for a,b in indices],
                                 [h[a][b].derivative(0).derivative(0).value for a,b in indices]])
                sums[side]+=weights[:,node,None,:]*values[None,:,:]
        limits=[]
        for position,fraction in enumerate((1.,.1,0.)):
            hv=[]
            dv=[]
            for side,sign in enumerate((-1.,1.)):
                d=sign*args.epsilon*(fraction-1.)
                hv.append(sums[side,:,0]+d*sums[side,:,1]+d*d*sums[side,:,2]/2)
                dv.append(sums[side,:,1]+d*sums[side,:,2])
            value_jump=hv[1]-hv[0]
            derivative_jump=dv[1]-dv[0]
            if baseline is not None:
                value_jump+=decode(baseline[position]['value_jump'])
                derivative_jump+=decode(baseline[position]['derivative_jump'])
            limits.append(dict(separation=2*args.epsilon*fraction,
                               value_jump=encode(value_jump),
                               derivative_jump=encode(derivative_jump),
                               derivative_error=encode(derivative_jump-targets),
                               maximum_value_jump=float(np.max(abs(value_jump))),
                               maximum_derivative_error=float(np.max(abs(derivative_jump-targets)))))
        case=dict(ellmax=ell,limits=limits)
        cases.append(case)
        print(dict(ellmax=ell,maximum_limit_value_jump=limits[-1]['maximum_value_jump'],
                   maximum_limit_derivative_error=limits[-1]['maximum_derivative_error']),flush=True)
        temporary=out.with_suffix('.tmp')
        temporary.write_text(json.dumps(dict(status='all_component_projected_diagnostic',jet_order=8,smooth_polar_tests=args.smooth,r0=6.,a=.6,m=2,
                 quadrature=args.quadrature,epsilon=args.epsilon,components=indices,tests=['P22','P32'],
                 expected_derivative_jump=targets.tolist(),cases=cases),indent=2)+'\n')
        temporary.replace(out)


if __name__=='__main__':
    main()
