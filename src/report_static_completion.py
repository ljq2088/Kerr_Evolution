"""Record vacuum/trace/gauge checks for the unmatched static completion basis."""
import json
from pathlib import Path
import numpy as np
from lorenz_completion import completion_metric
from lorenz_tensor import trace,lorenz_constraint,linearized_einstein


def main():
    rows=[]
    for a in (0.,.6,.8771530275949366):
        for r in (3.5,9.):
            for mode in 'ABCDEFG':
                g,h,expected=completion_metric(r,1.1,a=a,mode=mode)
                gauge=np.array([v.value for v in lorenz_constraint(g,h)])
                einstein=np.array([[v.value for v in row] for row in linearized_einstein(g,h)])
                rows.append(dict(a=a,r=r,mode=mode,trace_error=float(abs(trace(g,h).value-expected.value)),
                    lorenz_error=float(np.max(abs(gauge))),einstein_error=float(np.max(abs(einstein)))))
    maxima={key:max(row[key] for row in rows) for key in ('trace_error','lorenz_error','einstein_error')}
    if max(maxima.values())>2e-8:
        raise RuntimeError(f'Completion identity failure: {maxima}')
    path=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/static_completion_validation.json'
    path.write_text(json.dumps(dict(status='vacuum_basis_only_not_physical_completion',
        reference='https://arxiv.org/abs/2306.16459v3',y_radial_data='zero at r=6; homogeneous gauge freedom retained',
        missing=['physical boundary adaptation','particle matching coefficients','independent conserved-charge check'],
        maxima=maxima,rows=rows),indent=2)+'\n')
    print(json.dumps(maxima))


if __name__=='__main__':
    main()
