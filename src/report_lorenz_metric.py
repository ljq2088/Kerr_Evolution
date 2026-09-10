"""Audit nonstatic metric assembly; report failed gates without disguising them."""
import json
from pathlib import Path
import numpy as np
from lorenz_metric import nonstatic_metric
from lorenz_tensor import trace,lorenz_constraint,linearized_einstein
from lorenz_kappa import trace_field_jet


def main():
    cases=[]
    for r in (4.5,8.):
        g,h=nonstatic_metric(r,1.1,6.)
        _,hm=nonstatic_metric(r,1.1,6.,m=-2)
        values=np.array([[v.value for v in row] for row in h])
        minus=np.array([[v.value for v in row] for row in hm])
        scale=float(np.max(abs(values)))
        expected=trace_field_jet(g,6.,2).value
        case=dict(r=r,a=.6,r0=6.,ell=2,m=2,theta=1.1,
                  maximum_component=scale,
                  trace_error=float(abs(trace(g,h).value-expected)),
                  relative_lorenz=float(max(abs(v.value) for v in lorenz_constraint(g,h))/scale),
                  relative_einstein=float(max(abs(v.value) for row in linearized_einstein(g,h) for v in row)/scale),
                  raw_complex_conjugacy_defect=float(np.max(abs(values-minus.conjugate()))/scale))
        cases.append(case)
        print(json.dumps(case),flush=True)
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/metric_development.json'
    out.write_text(json.dumps(dict(status='nonstatic_development_source_matching_not_validated',spin1_normalization='factor 2 verified by Maxwell circularity',conjugacy_note='Complex AAB intermediate; conjugacy defect is not a physical-reality acceptance gate.',cases=cases),indent=2)+'\n')


if __name__=='__main__':
    main()
