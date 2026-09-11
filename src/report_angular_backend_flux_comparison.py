"""Compare complete source and flux channels from two angular eigensolvers."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('coverage',type=Path)
    parser.add_argument('dense_batch',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    coverage=json.loads(args.coverage.read_text())
    batch=json.loads(args.dense_batch.read_text())
    if batch['status']!='batch_completed_finite_resolution_not_converged':
        raise ValueError('Dense batch must be completed')
    available={(r['ell'],r['m']):r for r in coverage['infinity']['modes']}
    rows=[]
    for channel in batch['channels']:
        key=(channel['scalar_ell'],channel['scalar_m'])
        reference=available[key]
        if 'flux' not in reference:raise ValueError('Missing default channel')
        paths=[args.coverage.parent/reference['file'],args.dense_batch.parent/channel['file']]
        data=[json.loads(p.read_text()) for p in paths]
        for d in data:
            if d['status']!='truncated_single_mode_not_converged' or not d['propagating']:
                raise ValueError('Require completed propagating channels')
            if d.get('flux_validity')=='historical_unreliable_boundary_result':
                raise ValueError('Unreliable historical boundary')
        parameters=[copy.deepcopy(d['parameters']) for d in data]
        backends=[p['metric'].pop('angular_backend',None) for p in parameters]
        if backends!=[None,'dense-real-evd-diagnostic'] or parameters[0]!=parameters[1]:
            raise ValueError('Only the angular backend may differ')
        grids=[[(s['r'],s['weight']) for s in d['samples']] for d in data]
        if grids[0]!=grids[1]:raise ValueError('Different source nodes or weights')
        old,new=[np.array([complex(*s['source']) for s in d['samples']]) for d in data]
        weights=np.array([s['weight'] for s in data[0]['samples']])
        if np.any(weights<=0):raise ValueError('Require positive quadrature weights')
        fluxes={}
        for boundary,zkey in (('infinity','z_inf'),('horizon','z_h')):
            before,after=[d['flux'][boundary]['orbital_energy'] for d in data]
            z0,z1=[complex(*d[zkey]) for d in data]
            relative=after/before-1 if before else None
            if before:
                np.testing.assert_allclose(relative,abs(z1/z0)**2-1,atol=2e-14,rtol=1e-10)
            fluxes[boundary]=dict(default=before,dense=after,relative_change=relative)
        rows.append(dict(ell=key[0],m=key[1],
            inputs={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
            source_nodes=len(old),
            source_max_difference_over_max_source=float(max(abs(new-old))/max(abs(old))),
            source_weighted_relative_l2=float(np.sqrt(np.dot(weights,abs(new-old)**2)/np.dot(weights,abs(old)**2))),
            fluxes=fluxes))
    result=dict(status='complete_selected_channels_backend_audit_not_global_convergence',
        rows=rows,limitations=['Same finite source grid and metric cutoff',
            'Weighted source norm is a discrete diagnostic, not a continuum error bound',
            'Four channels do not establish accuracy of all modes or resolve paper normalization'])
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(rows),flush=True)


if __name__=='__main__':main()
