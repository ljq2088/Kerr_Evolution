"""Propagate only matching degrees outside the paper's cutoff guard.

Reuse versioned legacy jumps already recorded by the all-degree diagnostic;
never modify its cache or extend it while comparing cutoff resolutions.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from environment_source import ThresholdCloud
from environment_angular_diagnostic import install_dense_angular_diagnostic
from lorenz_jet import Jet
from report_full_gauge_closure import enc,decode,corrections
from source_provenance import local_dependency_hashes


def main():
    p=argparse.ArgumentParser();p.add_argument('--ellmax',type=int,default=10);p.add_argument('--quadrature',type=int,default=32)
    args=p.parse_args();root=Path(__file__).resolve().parents[1];folder=root/'docs/environment_reproduction'
    paths=[folder/f'all_component_matching_L{args.ellmax}_q{args.quadrature}_precise.json',
           folder/'full_gauge_closure_L8_precise.json',folder/'fresh_20260915_L18_mg-1_sl0.json']
    match,old,baseline=[json.loads(path.read_text()) for path in paths]
    for name,digest in match['implementation_sha256'].items():
        if hashlib.sha256((root/'src'/name).read_bytes()).hexdigest()!=digest:raise ValueError('Matching implementation differs: '+name)
    keep=args.ellmax-7
    if not 1<=keep<=len(old['rows']):raise ValueError('Need validated degrees and recorded legacy jumps')
    current=decode(match['solution']).reshape(args.ellmax,4)[:keep]
    legacy=np.array([decode(row['old_jumps']) for row in old['rows'][:keep]])
    install_dense_angular_diagnostic();Jet.coefficient_dtype=np.clongdouble
    cloud=ThresholdCloud();params=baseline['parameters'];Z=complex(*baseline['z_h'])
    dz=corrections(cloud,20.,current-legacy,params['source_panels'][0],params['source_panels'][-1]);total=np.sum(dz)
    result=dict(status='cutoff_guarded_gauge_correction_only_not_total_flux_validation',
        matching_L=args.ellmax,retained_ells=list(range(1,keep+1)),guard='ellmax - ell > 6',
        input_sha256={str(path.relative_to(root)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
        implementation_sha256=local_dependency_hashes(root/'src',['report_trusted_gauge_closure']),
        baseline_z_h=enc(Z),delta_z_h=enc(total),per_ell_delta_z_h=enc(dz),
        relative_flux_change=float(abs(1+total/Z)**2-1),
        limitations=['Only cutoff-guarded homogeneous gauge differences are changed',
          'Shared physical curvature/trace inputs and finite radial endpoints',
          'No replacement of all L18 modes or a new summed Fig.2 flux'])
    out=folder/f'trusted_gauge_closure_L{args.ellmax}.json';out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
