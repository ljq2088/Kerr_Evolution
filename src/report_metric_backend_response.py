"""Rebuild actual finite-grid scalar source with explicit Teukolsky methods.

A separate process is required per policy. No old tensor or scalar-source
cache is used. This distinguishes normalized radial-function comparisons from
the source, reconstruction and scalar-response combination actually used.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from environment_radial_backend_diagnostic import install,CALLS
from source_provenance import local_dependency_hashes


def main():
    p=argparse.ArgumentParser();p.add_argument('--spin2',default='AUTO');p.add_argument('--gauge',default='AUTO')
    p.add_argument('--rtol',type=float);p.add_argument('--ellmax',type=int,default=1)
    args=p.parse_args();install(args.spin2,args.gauge,args.rtol)
    from environment_source import ThresholdCloud,project_source
    from environment_lorenz_mode import LorenzMetricMode,ConjugateMetricMode
    from environment_radial import RadialGreen
    from environment_cloud import mode_flux
    from environment_angular_diagnostic import install_dense_angular_diagnostic
    # Deterministic angular input is fixed across all compared radial policies.
    install_dense_angular_diagnostic()
    root=Path(__file__).resolve().parents[1];folder=root/'docs/environment_reproduction'
    reference=folder/'fresh_20260915_L18_mg-1_sl0.json';ref=json.loads(reference.read_text())
    grid=ref['samples'];radii=np.array([s['r'] for s in grid]);weights=np.array([s['weight'] for s in grid])
    cloud=ThresholdCloud();metric=ConjugateMetricMode(LorenzMetricMode(20.,cloud.a,1,args.ellmax))
    green=RadialGreen(cloud.a,cloud.mu,cloud.omega+metric.omega,0,0,rmax=1000.,offset=1e-4,rtol=1e-11)
    encode=lambda z:[float(z.real),float(z.imag)]
    suffix=f'L{args.ellmax}_s2{args.spin2}_g{args.gauge}_rtol{args.rtol}'
    out=folder/f'metric_backend_response_{suffix}.json'
    if out.exists():raise FileExistsError('Choose a new policy/run file; this diagnostic must sample fresh')
    import pybhpt.radial,pybhpt.teuk,cybhpt_full
    artifacts={mod.__name__:dict(path=mod.__file__,sha256=hashlib.sha256(Path(mod.__file__).read_bytes()).hexdigest()) for mod in (pybhpt.radial,pybhpt.teuk,cybhpt_full)}
    start=time.perf_counter();samples=[]
    d=dict(status='sampling',ellmax=args.ellmax,spin2_method=args.spin2,gauge_method=args.gauge,rtol=args.rtol,
        angular_backend='DenseRealHarmonic',grid_reference_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),
        package_artifacts=artifacts,implementation_sha256=local_dependency_hashes(root/'src',['report_metric_backend_response']),samples=samples,
        limitations=['Finite 88x18 grid and specified metric ell cutoff; not a new full L18 total',
            'Trace/kappa resolvent and massive scalar Green solver unchanged',
            'Weyl source interface accepts method but does not expose radial ODE rtol'])
    def save():
        temp=out.with_suffix('.tmp');temp.write_text(json.dumps(d,indent=2)+'\n');temp.replace(out)
    try:
        for i,r in enumerate(radii):
            omega,source=project_source(cloud,[r],20.,0,0,metric,ntheta=18)
            if not np.all(np.isfinite(source)):raise ArithmeticError('Nonfinite projected source')
            samples.append(dict(r=float(r),weight=float(weights[i]),source=encode(source[0])))
            save();print(f'{i+1}/{len(radii)} r={r:.8g}',flush=True)
        J=np.array([complex(*s['source']) for s in samples]);Z=np.dot(weights,green.upsol.sol(radii)[0]*J)/green.w0
        d.update(status='finite_grid_backend_response_complete',z_h=encode(Z),flux=mode_flux(omega,0,cloud.omega,cloud.m,cloud.mu,cloud.a,0j,Z),
            elapsed_seconds=time.perf_counter()-start,radial_calls=CALLS)
        print('Z_H',Z,flush=True)
    except Exception as error:
        d.update(status='backend_or_sampling_failed',error=repr(error),radial_calls=CALLS,elapsed_seconds=time.perf_counter()-start);save();raise
    save()


if __name__=='__main__':main()
