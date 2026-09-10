"""Finite angular-cutoff continuity diagnostic away from the particle.

Neither a single cutoff nor nonzero radial separation establishes a source
jump condition. Save the sequence to guide the next convergence check.
"""
import json
from pathlib import Path
import numpy as np
from lorenz_metric import nonstatic_metric


def main():
    sums=[np.zeros((4,4),complex),np.zeros((4,4),complex)]
    cases=[]
    epsilon=5e-4
    for ell in range(2,7):
        for side,r in enumerate((6.-epsilon,6.+epsilon)):
            _,hp=nonstatic_metric(r,1.1,6.,ell=ell,order=6)
            _,hm=nonstatic_metric(r,1.1,6.,ell=ell,m=-2,order=6)
            p=np.array([[v.value for v in row] for row in hp])
            n=np.array([[v.value for v in row] for row in hm])
            sums[side]+=(p+n.conjugate())/2
        scale=max(float(np.max(abs(x))) for x in sums)
        case=dict(ellmax=ell,relative_two_sided_difference=float(np.max(abs(sums[1]-sums[0]))/scale),
                  maximum_component=scale)
        cases.append(case)
        print(case,flush=True)
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/metric_matching_development.json'
    out.write_text(json.dumps(dict(status='finite_cutoff_diagnostic_not_source_validation',r0=6.,a=.6,m=2,theta=1.1,radial_separation=2*epsilon,cases=cases),indent=2)+'\n')


if __name__=='__main__':
    main()
