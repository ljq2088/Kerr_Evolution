"""Compare auxiliary resolvent boundaries on sourced odd-ell m=1 modes."""
import json
import numpy as np
from pathlib import Path
from lorenz_kappa import kappa_jet
from lorenz_ghp import KerrGHP


def main():
    rows=[]
    a=.8771530275949366
    r0=41.6
    omega=1/(r0**1.5+a)
    output=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/kappa_boundary_comparison.json'
    # Even ell has zero equatorial scalar trace source at m=1.
    for ell in (1,7,9):
        for r in (20.,60.):
            g=KerrGHP(r,1.1,a,omega=omega,m=1,order=4)
            for method,rmax in [('series',2000.),('series',32000.),('coulomb',2000.),('coulomb',8000.),('coulomb',32000.)]:
                k=kappa_jet(g,r0,ell,rmax=rmax,infinity_method=method)
                values=[k.value,k.derivative(0).value,k.derivative(0).derivative(0).value]
                rows.append(dict(ell=ell,r=r,method=method,rmax=rmax,
                                 radial_jet=[[v.real,v.imag] for v in values]))
            print(f'ell={ell}, r={r}: completed',flush=True)
            output.write_text(json.dumps(dict(status='auxiliary_boundary_diagnostic_not_source_validation',
                parameters=dict(a=a,r0=r0,m=1,theta=1.1),rows=rows),indent=2)+'\n')
    summary=[]
    for ell in (1,7,9):
        for r in (20.,60.):
            group=[row for row in rows if row['ell']==ell and row['r']==r]
            values=lambda row:np.array([complex(*v) for v in row['radial_jet']])
            reference=values(group[-1])
            scale=float(np.max(abs(reference)))
            entry=dict(ell=ell,r=r,reference_scale=scale)
            for row,label in zip(group[:-1],('series2000','series32000','coulomb2000','coulomb8000')):
                entry[label+'_relative']=float(np.max(abs(values(row)-reference))/scale)
            summary.append(entry)
    output.write_text(json.dumps(dict(status='auxiliary_boundary_diagnostic_not_source_validation',
        parameters=dict(a=a,r0=r0,m=1,theta=1.1),rows=rows,summary=summary),indent=2)+'\n')


if __name__=='__main__':
    main()
