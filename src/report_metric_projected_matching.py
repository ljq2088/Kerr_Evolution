"""Angularly projected metric continuity, retaining per-spin contributions."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import lpmv
from lorenz_metric import spin2_metric,spin1_metric,spin0_metric


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--quadrature',type=int,default=8)
    parser.add_argument('--ellmax',type=int,default=6)
    args=parser.parse_args()
    x,w=np.polynomial.legendre.leggauss(args.quadrature)
    # Smooth scalar test function proportional to Y_22(theta,phi=0).
    weights=w*lpmv(2,2,x)
    sums=np.zeros((2,3),complex)
    cases=[]
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'/f'metric_projected_matching_q{args.quadrature}.json'
    for ell in range(2,args.ellmax+1):
        pieces=np.zeros((2,3),complex)
        for side,r in enumerate((6.-5e-4,6.+5e-4)):
            for theta,weight in zip(np.arccos(x),weights):
                for part,fun in enumerate((spin2_metric,spin1_metric,spin0_metric)):
                    _,hp=fun(r,theta,6.,ell=ell,order=6)
                    _,hm=fun(r,theta,6.,ell=ell,m=-2,order=6)
                    pieces[side,part]+=weight*(hp[0][0].value+hm[0][0].value.conjugate())/2
        sums+=pieces
        def encode(z):
            return [float(z.real),float(z.imag)]
        case=dict(ellmax=ell,tt_inside=encode(sum(sums[0])),tt_outside=encode(sum(sums[1])),
                  cumulative_piece_differences=[encode(v) for v in sums[1]-sums[0]],
                  relative_difference=float(abs(sum(sums[1]-sums[0]))/max(abs(sum(sums[0])),abs(sum(sums[1])))))
        cases.append(case)
        print(case,flush=True)
        out.write_text(json.dumps(dict(status='projected_diagnostic_not_full_source_validation',quadrature=args.quadrature,
                          r0=6.,a=.6,m=2,radial_separation=.001,test_function='P_2^2(cos theta)',cases=cases),indent=2)+'\n')


if __name__=='__main__':
    main()
