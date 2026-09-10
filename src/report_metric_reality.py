"""Complex AAB conjugacy diagnostic, not a physical-reality acceptance test.

The paper explicitly takes the real part of the complex reconstruction.
Raw +/-m conjugacy need not hold, even after angular summation.
"""
import json
from pathlib import Path
import numpy as np
from lorenz_metric import nonstatic_metric


def main():
    plus=np.zeros((4,4),complex)
    minus=plus.copy()
    cases=[]
    for ell in range(2,7):
        _,h=nonstatic_metric(8.,1.1,6.,ell=ell,order=8)
        _,hm=nonstatic_metric(8.,1.1,6.,ell=ell,m=-2,order=8)
        plus+=np.array([[v.value for v in row] for row in h])
        minus+=np.array([[v.value for v in row] for row in hm])
        error=float(np.max(abs(plus-minus.conjugate()))/np.max(abs(plus)))
        cases.append(dict(ellmax=ell,raw_complex_conjugacy_defect=error,maximum_component=float(np.max(abs(plus)))))
        print(cases[-1],flush=True)
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/metric_reality_development.json'
    out.write_text(json.dumps(dict(status='diagnostic_only_not_a_reality_gate',physical_projection='(h_m+conj(h_minus_m))/2; source matching still required',spin1_normalization='factor 2 verified by Maxwell circularity',r=8.,theta=1.1,r0=6.,a=.6,m=2,cases=cases),indent=2)+'\n')


if __name__=='__main__':
    main()
