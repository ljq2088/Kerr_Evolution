"""Angularly projected metric continuity, retaining per-spin contributions."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import lpmv
from lorenz_metric import spin2_metric,spin1_metric,spin0_metric
from environment_source import kerr_metric


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--quadrature',type=int,default=8)
    parser.add_argument('--ellmax',type=int,default=6)
    parser.add_argument('--epsilon',type=float,default=5e-4)
    args=parser.parse_args()
    metric=kerr_metric(6.,np.pi/2,.6)
    op=1/(6**1.5+.6)
    ut=1/np.sqrt(-metric[0,0]-2*op*metric[0,3]-op*op*metric[3,3])
    u_t=ut*(metric[0,0]+op*metric[0,3])
    expected_jump=-24*(u_t*u_t+metric[0,0]/2)/(ut*(36-12+.36))
    x,w=np.polynomial.legendre.leggauss(args.quadrature)
    # Smooth scalar test function proportional to Y_22(theta,phi=0).
    weights=w*lpmv(2,2,x)
    sums=np.zeros((2,3),complex)
    derivative_sums=sums.copy()
    cases=[]
    suffix='' if args.epsilon==5e-4 else f'_eps{args.epsilon:g}'
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'/f'metric_chiral_matching_q{args.quadrature}{suffix}.json'
    for ell in range(2,args.ellmax+1):
        pieces=np.zeros((2,3),complex)
        derivatives=pieces.copy()
        for side,r in enumerate((6.-args.epsilon,6.+args.epsilon)):
            for theta,weight in zip(np.arccos(x),weights):
                for part,fun in enumerate((spin2_metric,spin1_metric,spin0_metric)):
                    extra={'full_current':True} if part==1 else {}
                    _,hp=fun(r,theta,6.,ell=ell,order=6,**extra)
                    _,hm=fun(r,theta,6.,ell=ell,m=-2,order=6,**extra)
                    pieces[side,part]+=weight*(hp[0][0].value+hm[0][0].value.conjugate())/2
                    derivatives[side,part]+=weight*(hp[0][0].derivative(0).value+hm[0][0].derivative(0).value.conjugate())/2
        sums+=pieces
        derivative_sums+=derivatives
        def encode(z):
            return [float(z.real),float(z.imag)]
        case=dict(ellmax=ell,tt_inside=encode(sum(sums[0])),tt_outside=encode(sum(sums[1])),
                  cumulative_piece_differences=[encode(v) for v in sums[1]-sums[0]],
                  derivative_jump=encode(sum(derivative_sums[1]-derivative_sums[0])),
                  relative_difference=float(abs(sum(sums[1]-sums[0]))/max(abs(sum(sums[0])),abs(sum(sums[1])))))
        cases.append(case)
        print(case,flush=True)
        out.write_text(json.dumps(dict(status='projected_diagnostic_not_full_source_validation',quadrature=args.quadrature,
                          r0=6.,a=.6,m=2,radial_separation=2*args.epsilon,expected_derivative_jump=float(expected_jump),test_function='P_2^2(cos theta)',cases=cases),indent=2)+'\n')


if __name__=='__main__':
    main()
