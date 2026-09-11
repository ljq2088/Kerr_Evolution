"""Direct circular-orbit gravitational flux sums with both signs of m."""
import argparse
import importlib.metadata
import json
from pathlib import Path
import numpy as np
from pybhpt.geo import KerrGeodesic
from pybhpt.teuk import TeukolskyMode
from pybhpt.flux import FluxMode


def mode_flux(a,r,ell,m):
    geo=KerrGeodesic(a,r,0.,1.)
    mode=TeukolskyMode(-2,ell,m,0,0,geo)
    mode.solve(geo)
    flux=FluxMode(geo,mode)
    omega_p=1/(r**1.5+a)
    for boundary in ('I','H'):
        np.testing.assert_allclose(flux.energy[boundary],omega_p*flux.angularmomentum[boundary],
                                   rtol=1e-12,atol=1e-300)
    return dict(infinity=flux.energy['I'],horizon=flux.energy['H'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--a',type=float,default=.8771530275949366)
    parser.add_argument('--radius',type=float,default=20.)
    parser.add_argument('--ellmax',type=int,default=12)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.ellmax<2 or not 0<=args.a<1:raise ValueError('Invalid mode range or spin')
    benchmark=mode_flux(.9,10.,2,2)
    expected=dict(infinity=.000022273,horizon=-5.9836e-8)
    for key in expected:
        np.testing.assert_allclose(benchmark[key],expected[key],rtol=3e-5)
    rows=[];shells=[];totals=dict(infinity=0.,horizon=0.)
    result=dict(status='in_progress',parameters=dict(a=args.a,rp=args.radius,ellmax=args.ellmax,
        spinweight=-2,k=0,n=0,eccentricity=0.,inclination_x=1.,units='per q^2'),
        package_version=importlib.metadata.version('pybhpt'),
        benchmark=dict(source='https://bhptoolkit.org/modules/teukolsky/',
            parameters=dict(a=.9,rp=10.,ell=2,m=2),rounded_reference=expected,computed=benchmark),
        modes=rows,shells=shells,
        limitations=['Both nonzero m signs explicitly solved; static m=0 carries no radiation',
            'Finite ell sum; shell changes are not a rigorous tail bound',
            'Independent public Toolkit rounded benchmark checks conventions, not full accuracy'])
    for ell in range(2,args.ellmax+1):
        shell=dict(infinity=0.,horizon=0.)
        for m in range(-ell,ell+1):
            if m==0:continue
            flux=mode_flux(args.a,args.radius,ell,m)
            if not all(np.isfinite(v) for v in flux.values()):raise ValueError('Nonfinite flux')
            rows.append(dict(ell=ell,m=m,**flux))
            for key in shell:shell[key]+=flux[key]
        for key in shell:totals[key]+=shell[key]
        shells.append(dict(ell=ell,contribution=shell,cumulative=dict(totals)))
        print(shells[-1],flush=True)
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    result.update(status='finite_gravitational_flux_sum_not_full_paper_reproduction',totals=totals)
    args.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
