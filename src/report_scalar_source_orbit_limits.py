"""Direct one-sided scalar-source projections at the particle orbit.

Unlike tensor matching with selected polynomial tests, weights here are the
actual cloud Hessian and scalar spheroidal harmonics. Finite-L discontinuities
are diagnostic, not a proof of a distributional source error.
"""
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np


def evaluate(offset):
    from environment_source import ThresholdCloud, angular_mode
    from environment_angular_diagnostic import install_dense_angular_diagnostic
    from lorenz_metric import nonstatic_metric
    install_dense_angular_diagnostic()
    cloud=ThresholdCloud(alpha=.3);r0=20.;r=r0+offset
    omega=cloud.omega-1/(r0**1.5+cloud.a)
    x,w=np.polynomial.legendre.leggauss(18);theta=np.arccos(x)
    ells=(0,1,2)
    tests=np.array([angular_mode(theta,l,0,cloud.a**2*(omega**2-cloud.mu**2))[0] for l in ells])
    weights=2*np.pi*w*(r*r+cloud.a**2*x*x)*tests
    kernels=[]
    for t in theta:
        _,hessian,inv=cloud.hessian(r,t);kernels.append(inv@hessian@inv)
    cumulative=np.zeros(3,complex);rows=[]
    folder=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    path=folder/f'scalar_source_orbit_limit_offset{offset:+g}.json'
    def enc(a):return np.stack([a.real,a.imag],axis=-1).tolist()
    for L in range(1,19):
        samples=[]
        for t,kernel in zip(theta,kernels):
            _,h=nonstatic_metric(r,float(t),r0,cloud.a,L,1,order=6)
            values=np.array([[v.value for v in row] for row in h]).conjugate()
            samples.append(np.einsum('ab,ab->',values,kernel))
        piece=weights@np.array(samples);cumulative+=piece
        rows.append(dict(metric_ellmax=L,piece=enc(piece),cumulative=enc(cumulative)))
        result=dict(status='one_sided_source_in_progress' if L<18 else 'one_sided_source_completed',
            r0=r0,offset=offset,a=cloud.a,alpha=.3,scalar_m=0,scalar_ells=ells,
            quadrature=18,angular_backend='dense-real',kappa_backend='production-finite-difference',
            metric_jet_order=6,rows=rows)
        tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(result,indent=2)+'\n');tmp.replace(path)
        print('offset',offset,'metric ellmax',L,flush=True)
    return path.name


def main():
    with ProcessPoolExecutor(max_workers=4) as pool:
        jobs=[pool.submit(evaluate,offset) for offset in (-.001,-.0005,.0005,.001)]
        for job in as_completed(jobs):print('Completed',job.result(),flush=True)


if __name__=='__main__':main()
