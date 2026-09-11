"""Discrete integral conditioning and tighter homogeneous radial integration.

No source samples, source cutoffs, or angular factors are modified.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from environment_radial import RadialGreen


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs',type=Path,nargs='+')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();rows=[]
    for path in args.inputs:
        data=json.loads(path.read_text());p=data['parameters'];metric=p['metric']
        if data['status']!='truncated_single_mode_not_converged' or not data['propagating']:
            raise ValueError('Require completed propagating source')
        if data.get('flux_validity')=='historical_unreliable_boundary_result':raise ValueError('Invalid old boundary')
        r=np.array([s['r'] for s in data['samples']]);w=np.array([s['weight'] for s in data['samples']])
        source=np.array([complex(*s['source']) for s in data['samples']]);runs=[]
        for tolerance in (1e-11,1e-13):
            green=RadialGreen(metric['a'],p['alpha'],p['omega'],p['scalar_ell'],p['scalar_m'],
                rmax=p['green_outer_radius'],offset=p['green_horizon_offset'],rtol=tolerance,
                infinity_method=p.get('infinity_method','series'))
            terms=w*green.insol.sol(r)[0]*source/green.w0
            amplitude=terms.sum()
            runs.append(dict(rtol=tolerance,amplitude=[amplitude.real,amplitude.imag],
                discrete_condition_number=float(sum(abs(terms))/abs(amplitude)),
                wronskian_relative_spread=float(max(abs(green.wronskian(r)/green.w0-1)))))
        baseline=complex(*runs[0]['amplitude']);tight=complex(*runs[1]['amplitude'])
        recorded=complex(*data['z_inf'])
        if abs(baseline-recorded)>1e-10*abs(recorded):raise ValueError('Baseline differs from recorded amplitude')
        row=dict(input=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            ell=p['scalar_ell'],m=p['scalar_m'],runs=runs,
            relative_complex_amplitude_change=float(abs(tight-baseline)/abs(baseline)),
            relative_infinity_flux_change=float(abs(tight/baseline)**2-1))
        rows.append(row);print(row,flush=True)
    result=dict(status='same_source_radial_integrator_audit_not_full_convergence',rows=rows,
        limitations=['Condition number belongs to the finite sampled sum, not a continuum error bound',
            'No source-grid, angular, metric, or finite-boundary refinement in this check',
            'Near-unit conditioning does not bound unknown absolute source errors'])
    args.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
