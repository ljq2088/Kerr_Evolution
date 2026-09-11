"""Actual Lorenz-source angular selection audit in the spherical limit."""
import json
from pathlib import Path
import numpy as np
from environment_schwarzschild_cloud import SchwarzschildCloud
from environment_source import angular_mode
from lorenz_metric import nonstatic_metric


def main():
    cloud=SchwarzschildCloud(alpha=.3,freeze_decay=True)
    r,r0,m=10.,20.,2
    x,w=np.polynomial.legendre.leggauss(12)
    theta=np.arccos(x)
    angular={ell:angular_mode(theta,ell,3,0.)[0] for ell in (3,4)}
    kernels=[]
    for t in theta:
        _,H,inverse=cloud.hessian(r,t)
        kernels.append(inverse@H@inverse)
    rows=[]
    for ell in range(2,7):
        values=[]
        for t,kernel in zip(theta,kernels):
            _,h=nonstatic_metric(r,t,r0,a=0.,ell=ell,m=m,order=6)
            values.append(np.einsum('ij,ij->',np.array([[v.value for v in row] for row in h]),kernel))
        projections={str(l):2*np.pi*r*r*np.dot(w,angular[l]*values) for l in angular}
        row=dict(metric_ell=ell,source={l:[v.real,v.imag] for l,v in projections.items()})
        rows.append(row)
        print(row,flush=True)
    reference=abs(sum(complex(*row['source']['3']) for row in rows if row['metric_ell']<=4))
    forbidden_parity=max(abs(complex(*row['source']['4'])) for row in rows)
    forbidden_triangle=max(abs(complex(*row['source']['3'])) for row in rows if row['metric_ell']>4)
    result=dict(status='single_radius_selection_diagnostic',parameters=dict(alpha=.3,r=r,r0=r0,metric_m=2,angular_order=12),
        source_scale=reference,opposite_parity_relative=forbidden_parity/reference,
        outside_triangle_relative=forbidden_triangle/reference,rows=rows)
    path=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/schwarzschild_source_selection.json'
    path.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    main()
