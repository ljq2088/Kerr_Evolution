"""Independent coordinate-curvature check of AAB circularity."""
import json
from pathlib import Path
from lorenz_metric import spin2_metric,homogeneous_field_jet
from lorenz_weyl import weyl_amplitudes
from lorenz_tensor import extreme_weyl


def main():
    cases=[]
    for r in (4.5,8.):
        g,h=spin2_metric(r,1.1,6.)
        computed=extreme_weyl(g,h)
        bc,index=('In',1) if r<6 else ('Up',0)
        amplitudes=weyl_amplitudes(6.,.6,2,2)
        expected=[homogeneous_field_jet(g,s,2,amplitudes[s][index],bc).value for s in (2,-2)]
        case=dict(r=r,relative_error=[float(abs(v/w-1)) for v,w in zip(computed,expected)],
                  ratio=[[float((v/w).real),float((v/w).imag)] for v,w in zip(computed,expected)])
        print(case,flush=True)
        cases.append(case)
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/metric_curvature_development.json'
    out.write_text(json.dumps(dict(status='independent_curvature_audit',cases=cases),indent=2)+'\n')


if __name__=='__main__':
    main()
