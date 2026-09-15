"""Fresh dipole scalar response with independently MST-anchored gauge radial data.

Keeps the trace/kappa resolvent separate and records the unchanged baseline
implementation hashes. Run each spin selection in its own process.
"""
import argparse
import hashlib
import importlib
import json
from pathlib import Path
import time
import numpy as np
from pybhpt.radial import RadialTeukolsky as NativeRadial
from environment_mst_anchored_radial import MSTAnchoredRadial
from source_provenance import local_dependency_hashes

ROOT = Path(__file__).resolve().parents[1]


def install(spins, anchor, up_anchor, rtol):
    import lorenz_metric
    modules = [importlib.import_module(name) for name in
               ('lorenz_metric', 'lorenz_spin1', 'lorenz_spin1_chiral', 'lorenz_chi')]
    calls = []
    class SelectedRadial:
        def __init__(self, s, ell, m, a, omega, radii):
            self.input = dict(spin=s, ell=ell, m=m, a=a, omega=omega)
            self.selected = s in spins
            self.impl = (MSTAnchoredRadial(s, ell, m, a, omega, radii, anchor=anchor, up_anchor=up_anchor, rtol=rtol)
                         if self.selected else NativeRadial(s, ell, m, a, omega, radii))
        def __getattr__(self, name):
            return getattr(self.impl, name)
        def __call__(self, *args, **kwargs):
            return self.impl(*args, **kwargs)
        def solve(self, method='AUTO', bc=None, rtol=None):
            self.impl.solve(method=method, bc=bc, rtol=rtol)
            for side in (('In', 'Up') if bc is None else (bc,)):
                values = np.r_[self.impl.radialsolutions(side), self.impl.radialderivatives(side)]
                if not np.all(np.isfinite(values)) or not np.any(values):
                    raise ArithmeticError(f'Invalid radial output {self.input}, {side}')
            calls.append(dict(self.input, backend='MST_anchored' if self.selected else 'AUTO',
                              bc=bc, diagnostics=getattr(self.impl, 'diagnostics', None)))
    for module in modules:
        for obj in vars(module).values():
            if callable(getattr(obj, 'cache_clear', None)):
                obj.cache_clear()
        module.RadialTeukolsky = SelectedRadial
    return calls


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--spins', choices=['all', 'zero', 'one'], default='all')
    parser.add_argument('--anchor', type=float, default=3.)
    parser.add_argument('--up-anchor', type=float, default=20.)
    parser.add_argument('--rtol', type=float, default=3e-14)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    out = Path(args.output)
    if out.exists():
        raise FileExistsError('A fresh output path is required')
    baseline_path = ROOT/'docs/environment_reproduction/metric_backend_response_L1_s2AUTO_gAUTO_rtolNone.json'
    baseline = json.loads(baseline_path.read_text())
    # Reusing comparison numbers is safe only with their exact original producer.
    for name, expected in baseline['implementation_sha256'].items():
        actual = hashlib.sha256((ROOT/'src'/f'{name}.py').read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f'Baseline implementation changed: {name}')
    for artifact in baseline['package_artifacts'].values():
        if hashlib.sha256(Path(artifact['path']).read_bytes()).hexdigest() != artifact['sha256']:
            raise ValueError('Baseline external radial library changed')
    spins = {'all': (0,-1,1), 'zero': (0,), 'one': (-1,1)}[args.spins]
    calls = install(spins, args.anchor, args.up_anchor, args.rtol)
    from environment_source import ThresholdCloud, project_source
    from environment_lorenz_mode import LorenzMetricMode, ConjugateMetricMode
    from environment_radial import RadialGreen
    from environment_cloud import mode_flux
    from environment_angular_diagnostic import install_dense_angular_diagnostic
    install_dense_angular_diagnostic()
    cloud = ThresholdCloud()
    metric = ConjugateMetricMode(LorenzMetricMode(20., cloud.a, 1, 1))
    green = RadialGreen(cloud.a,cloud.mu,cloud.omega+metric.omega,0,0,rmax=1000.,offset=1e-4,rtol=1e-11)
    radii = np.array([s['r'] for s in baseline['samples']])
    weights = np.array([s['weight'] for s in baseline['samples']])
    encode = lambda z: [float(z.real), float(z.imag)]
    started = time.perf_counter()
    report = dict(status='sampling', spins=list(spins), anchor=args.anchor, up_anchor=args.up_anchor, rtol=args.rtol,
                  baseline_sha256=hashlib.sha256(baseline_path.read_bytes()).hexdigest(),
                  implementation_sha256=local_dependency_hashes(ROOT/'src', ['report_gauge_mst_source_response']),
                  baseline_implementation_verified=True, samples=[],
                  limitations=['Dipole metric ell=1 and scalar00, finite 88x18 grid',
                               'Trace/kappa, angular data and massive scalar Green are unchanged',
                               'MST seed replaces boundary initial data as well as radial integration',
                               'No production backend is changed'])
    def save():
        tmp=out.with_suffix('.tmp')
        tmp.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
        tmp.replace(out)
    try:
        for i, r in enumerate(radii):
            omega, source=project_source(cloud,[r],20.,0,0,metric,ntheta=18)
            if not np.all(np.isfinite(source)):
                raise ArithmeticError('Nonfinite projected scalar source')
            report['samples'].append(dict(r=float(r),weight=float(weights[i]),source=encode(source[0])))
            save()
            print(f'{i+1}/{len(radii)} r={r:.8g}', flush=True)
        source=np.array([complex(*s['source']) for s in report['samples']])
        previous=np.array([complex(*s['source']) for s in baseline['samples']])
        zh=np.dot(weights, green.upsol.sol(radii)[0]*source)/green.w0
        flux=mode_flux(omega,0,cloud.omega,cloud.m,cloud.mu,cloud.a,0j,zh)
        report.update(status='completed_mst_anchored_dipole_response', z_h=encode(zh), flux=flux,
                      z_h_relative_difference=abs(zh-complex(*baseline['z_h']))/abs(complex(*baseline['z_h'])),
                      horizon_flux_relative_difference=flux['horizon']['orbital_energy']/baseline['flux']['horizon']['orbital_energy']-1,
                      source_weighted_relative_l2=float(np.linalg.norm(np.sqrt(abs(weights))*(source-previous))/np.linalg.norm(np.sqrt(abs(weights))*previous)),
                      elapsed_seconds=time.perf_counter()-started, radial_calls=calls)
        print('Horizon relative difference',report['horizon_flux_relative_difference'],flush=True)
    except Exception as error:
        report.update(status='failed', error=repr(error), elapsed_seconds=time.perf_counter()-started)
        save()
        raise
    save()


if __name__ == '__main__':
    main()
