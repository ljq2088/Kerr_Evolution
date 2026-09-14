"""Check the scalar-potential map between the 2023 and 2024 conventions.

This is a vacuum identity, not a source-matching or full-gauge equivalence proof.
"""
import json
from pathlib import Path
import numpy as np
from lorenz_ghp import KerrGHP
from lorenz_spin1 import cky_tensor
from lorenz_tensor import vector_covariant_derivative
from lorenz_metric import homogeneous_field_jet
from lorenz_chi import chi_amplitudes
from environment_trace_variation import TraceMassVariation
from paper_jump_basis import scalar_components
from report_nonstatic_paper_projection import transform,PAIRS


def check_case(a,r0,m,ell,r):
    omega=m/(r0**1.5+a)
    g=KerrGHP(r,1.1,a,omega=omega,m=m,order=6)
    variation=TraceMassVariation(r0,a,ell,m)
    h,dh=variation.field_jet(g)
    k24=-1j*omega*dh
    bc,index=('In',1) if r<r0 else ('Up',0)
    chi=homogeneous_field_jet(g,0,ell,chi_amplitudes(r0,a,ell,m)[index],bc)
    k23=1j/(2*omega)*(k24-chi)
    f=cky_tensor(g)
    rawxi=[sum(f[i][j]*g.inv[j][k]*g.partial(h,k)/2 for j in range(4) for k in range(4))
           +g.partial(k24-chi,i) for i in range(4)]
    dxi=vector_covariant_derivative(g,rawxi)
    metric=[-1j/omega*(dxi[i][j]+dxi[j][i]) for i,j in PAIRS]
    T,dT=transform(r,1.1,a)
    vals=np.array([v.value for v in metric]);ders=np.array([v.derivative(0).value for v in metric])
    reference=(T@vals)[[0,2,4]];dreference=(T@ders+dT@vals)[[0,2,4]]
    direct=scalar_components(g,h,k23)
    err=np.max(abs(np.array([v.value for v in direct])-reference))/np.max(abs(reference))
    derr=np.max(abs(np.array([v.derivative(0).value for v in direct])-dreference))/np.max(abs(dreference))
    residual=abs(g.scalar_wave(k23).value-h.value/2)/max(abs(h.value)/2,1e-30)
    return dict(a=a,r0=r0,m=m,ell=ell,r=r,metric_relative_error=float(err),
                derivative_relative_error=float(derr),kappa_equation_relative_residual=float(residual))


def main():
    cases=[(a,r0,m,ell,r) for a,r0,m,ell,rr in
        ((0.,6.,1,1,(4.,9.)),(.8771530275949366,20.,1,1,(10.,40.)),(.6,6.,2,2,(4.,9.))) for r in rr]
    rows=[]
    for case in cases:
        row=check_case(*case);rows.append(row);print(json.dumps(row),flush=True)
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/paper_scalar_mapping_audit.json'
    out.write_text(json.dumps(dict(status='vacuum_scalar_sector_identity',
      map='kappa_2023 = i*(kappa_2024 - chi_compact)/(2*omega)',rows=rows,
      limitations=['Nonstatic vacuum points only', 'No point-source matching proof',
                   'No equivalence assertion for the full spin-2/spin-1 sector allocations',
                   'No change to the physical scalar flux']),indent=2)+'\n')


if __name__=='__main__':main()
