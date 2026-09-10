"""Recover Maxwell scalars from reconstructed vector, independently of gauge tests."""
import json
from pathlib import Path
from lorenz_metric import spin1_metric,homogeneous_field_jet
from lorenz_spin1 import spin1_amplitudes
from lorenz_tensor import vector_covariant_derivative


def main():
    cases=[]
    for r in (4.5,8.):
        g,xi=spin1_metric(r,1.1,6.,return_vector=True)
        d=vector_covariant_derivative(g,xi)
        F=[[d[j][i]-d[i][j] for j in range(4)] for i in range(4)]
        l,n,mm,mb=g.tetrad
        computed=[sum(F[i][j]*u[i]*v[j] for i in range(4) for j in range(4)).value
                  for u,v in ((l,mm),(mb,n))]
        bc,index=('In',1) if r<6 else ('Up',0)
        expected=[homogeneous_field_jet(g,s,2,spin1_amplitudes(6.,.6,2,2,s)[index],bc).value for s in (1,-1)]
        case=dict(r=r,ratio=[[float((v/w).real),float((v/w).imag)] for v,w in zip(computed,expected)])
        cases.append(case)
        print(case,flush=True)
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/spin1_circularity_development.json'
    out.write_text(json.dumps(dict(status='normalization_audit',cases=cases),indent=2)+'\n')


if __name__=='__main__':
    main()
