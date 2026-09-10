"""Separate angular quadrature and metric-ell truncation in actual scalar sources."""
import json
from pathlib import Path
import numpy as np
from environment_source import ThresholdCloud, angular_mode
from lorenz_metric import nonstatic_metric


def main():
    cloud=ThresholdCloud(alpha=.3)
    omega=cloud.omega+2/(20**1.5+cloud.a)
    path=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/source_truncation.json'
    parameters=dict(alpha=.3,rp=20.,scalar_ell=3,scalar_m=3,metric_m=2,
                    radii=[10.,30.],angular_orders=[6,10,14],metric_ellmax=8)
    report=dict(status='sampling',parameters=parameters,rows=[])
    if path.exists():
        report=json.loads(path.read_text())
        if report['parameters']!=parameters:
            raise ValueError('Existing audit parameters differ')
    done={(r['radius'],r['angular_order'],r['metric_ell']) for r in report['rows']}
    for r in parameters['radii']:
        for n in parameters['angular_orders']:
            x,w=np.polynomial.legendre.leggauss(n)
            theta=np.arccos(x)
            s=angular_mode(theta,3,3,cloud.a**2*(omega**2-cloud.mu**2))[0]
            kernels=[]
            for t in theta:
                _,H,inv=cloud.hessian(r,t)
                kernels.append(inv@H@inv)
            for ell in range(2,9):
                if (r,n,ell) in done:
                    continue
                values=[]
                for t,kernel in zip(theta,kernels):
                    _,h=nonstatic_metric(r,t,20.,cloud.a,ell,2,order=6)
                    h=np.array([[v.value for v in row] for row in h])
                    values.append(np.einsum('ij,ij->',h,kernel))
                value=2*np.pi*np.dot(w,s*(r*r+cloud.a**2*x*x)*values)
                report['rows'].append(dict(radius=r,angular_order=n,metric_ell=ell,
                                           source=[float(value.real),float(value.imag)]))
                temporary=path.with_suffix('.tmp')
                temporary.write_text(json.dumps(report,indent=2)+'\n')
                temporary.replace(path)
                print(f'r={r:g}, ntheta={n}, ell_g={ell}: {value}',flush=True)
    report['status']='pointwise_diagnostic_not_flux_convergence'
    summary=[]
    for r in parameters['radii']:
        sums={}
        for n in parameters['angular_orders']:
            pieces={row['metric_ell']:complex(*row['source']) for row in report['rows']
                    if row['radius']==r and row['angular_order']==n}
            for upper in (4,6,8):
                sums[n,upper]=sum(pieces[l] for l in range(2,upper+1))
        reference=sums[14,8]
        summary.append(dict(radius=r,
            angular_6_to_10_at_L8_relative=abs(sums[6,8]-sums[10,8])/abs(reference),
            angular_10_to_14_at_L8_relative=abs(sums[10,8]-reference)/abs(reference),
            metric_L4_to_L8_at_n14_relative=abs(sums[14,4]-reference)/abs(reference),
            metric_L6_to_L8_at_n14_relative=abs(sums[14,6]-reference)/abs(reference)))
    report['summary']=summary
    path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
