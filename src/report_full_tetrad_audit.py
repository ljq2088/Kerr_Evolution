"""Compare all ten independent explicit 2023 tetrad components."""
import json,hashlib
from pathlib import Path
import numpy as np
from paper_full_tetrad import PAIRS,scalar_tetrad,vector_tetrad,tensor_tetrad,project_metric,to_bl
from paper_jump_basis import gauge_basis,radial_jet,angular_jet
from lorenz_ghp import KerrGHP
from lorenz_spin1 import cky_tensor
from lorenz_tensor import vector_covariant_derivative
from lorenz_metric import spin2_metric,homogeneous_field_jet,_homogeneous_radial_data
from lorenz_weyl import weyl_amplitudes
from lorenz_chi import chi_amplitudes
from environment_trace_variation import TraceMassVariation
from environment_angular_diagnostic import install_dense_angular_diagnostic
from lorenz_jet import Jet


def compare(direct,reference,g):
    d=np.array([v.value for v in direct],complex);ref=project_metric(reference,float(g.r.value.real),float(g.theta.value.real),g.a)
    scale=max(np.max(abs(ref)),1e-300);norm=float(np.max(abs(d-ref))/scale)
    bl=to_bl(d,float(g.r.value.real),float(g.theta.value.real),g.a)
    component=np.abs(d-ref)/np.maximum(abs(ref),scale*1e-10)
    return dict(norm_relative_error=norm,component_relative_floor=component.tolist(),
        bl_relative_error=float(np.max(abs(bl-reference))/max(np.max(abs(reference)),1e-300)))


def main():
    install_dense_angular_diagnostic();Jet.coefficient_dtype=np.clongdouble;rows=[]
    for a,r0,m,ell in [(0.,6.,1,1),(.6,6.,2,2),(.8771530275949366,20.,1,1),(.8771530275949366,20.,1,2)]:
        variation=TraceMassVariation(r0,a,ell,m)
        for r in ([3.,10.,40.] if r0==20 else [4.,9.]):
            t=1.1;w=m/(r0**1.5+a);g=KerrGHP(r,t,a,omega=w,m=m,order=10)
            h,dh=variation.field_jet(g);bc,index=('In',1) if r<r0 else ('Up',0)
            chi=homogeneous_field_jet(g,0,ell,chi_amplitudes(r0,a,ell,m)[index],bc);k24=-1j*w*dh;k23=1j*(k24-chi)/(2*w)
            ky=cky_tensor(g);raw=[sum(ky[i][j]*g.inv[j][k]*g.partial(h,k)/2 for j in range(4) for k in range(4))+g.partial(k24-chi,i) for i in range(4)]
            dx=vector_covariant_derivative(g,raw);ref=np.array([[-1j/w*(dx[i][j].value+dx[j][i].value) for j in range(4)] for i in range(4)],complex)
            scalar=compare(scalar_tetrad(g,h,k23),ref,g)
            p0,p1=.31+.22j,-.11+.53j;basis=[gauge_basis(g,ell,'spin1',i,True) for i in (0,1)]
            vector=[p0*basis[0][i]+p1*basis[1][i] for i in range(4)]
            cov=[sum(g.g[i][j]*vector[j] for j in range(4)) for i in range(4)]
            dx=vector_covariant_derivative(g,cov);ref=np.array([[-dx[i][j].value-dx[j][i].value for j in range(4)] for i in range(4)],complex)
            spin1=compare(vector_tetrad(g,ell,p0,p1),ref,g);tensor=None
            if ell>=2:
                gt,ht=spin2_metric(r,t,r0,a,ell,m,order=10)
                _,R,Rp=_homogeneous_radial_data(-2,ell,m,a,w,r,bc);amp=-4*weyl_amplitudes(r0,a,ell,m)[-2][index]
                ref=np.array([[v.value for v in line] for line in ht],complex)
                tensor=compare(tensor_tetrad(gt,ell,amp*R,amp*Rp),ref,gt)
            row=dict(a=a,r0=r0,m=m,ell=ell,r=r,theta=t,scalar=scalar,spin1=spin1,spin2=tensor)
            rows.append(row);print(json.dumps(row),flush=True)
    root=Path(__file__).resolve().parents[1];out=root/'docs/environment_reproduction/full_tetrad_formula_audit.json'
    out.write_text(json.dumps(dict(status='ten_component_formula_check_shared_potential_input',rows=rows,
        implementation_sha256={n:hashlib.sha256((root/'src'/n).read_bytes()).hexdigest() for n in ['paper_full_tetrad.py','lorenz_metric.py']},
        limitations=['Independent explicit tetrad formulas; same curvature and trace inputs',
          'Spin1 uses arbitrary homogeneous potentials and an independent covariant gauge transformation',
          'Not a global sourced metric matching validation',
          'Extended jets, double radial/angular data; errors retain conditioning information']),indent=2)+'\n')


if __name__=='__main__':main()
