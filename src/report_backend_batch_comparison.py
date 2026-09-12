"""Compare complete default/dense batches including bound scalar channels."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import numpy as np


def compare(default_path,dense_path):
    batches=[json.loads(p.read_text()) for p in (default_path,dense_path)]
    if any(b['status']!='batch_completed_finite_resolution_not_converged' for b in batches):
        raise ValueError('Both batches must be complete')
    parameters=[copy.deepcopy(b['parameters']) for b in batches]
    if [p.pop('dense_angular',False) for p in parameters]!=[False,True] or parameters[0]!=parameters[1]:
        raise ValueError('Only the dense angular backend may differ')
    channels=[{(c['scalar_ell'],c['scalar_m']):c for c in b['channels']} for b in batches]
    if not channels[0] or channels[0].keys()!=channels[1].keys():
        raise ValueError('Nonempty matching channel sets required')
    rows=[]
    for mode in sorted(channels[0]):
        paths=[p.parent/c[mode]['file'] for p,c in zip((default_path,dense_path),channels)]
        data=[json.loads(p.read_text()) for p in paths]
        params=[copy.deepcopy(d['parameters']) for d in data]
        if [p['metric'].pop('angular_backend',None) for p in params]!=[None,'dense-real-evd-diagnostic'] or params[0]!=params[1]:
            raise ValueError('Response parameters differ beyond angular backend')
        if data[0]['propagating']!=data[1]['propagating']:
            raise ValueError('Propagation branches differ')
        for d,c in zip(data,channels):
            if d['status']!='truncated_single_mode_not_converged' or d['flux']!=c[mode]['flux']:
                raise ValueError('Invalid or inconsistent response')
            if d.get('flux_validity')=='historical_unreliable_boundary_result':
                raise ValueError('Unreliable historical boundary')
        grids=[[(s['r'],s['weight']) for s in d['samples']] for d in data]
        if grids[0]!=grids[1]:raise ValueError('Source quadrature differs')
        source=[np.array([complex(*s['source']) for s in d['samples']]) for d in data]
        weights=np.array([s['weight'] for s in data[0]['samples']])
        if not np.isfinite(source).all() or not np.isfinite(weights).all() or np.any(weights<=0):
            raise ValueError('Finite sources and positive finite weights required')
        flux={}
        for boundary,key in (('horizon','z_h'),('infinity','z_inf')):
            before,after=[d['flux'][boundary]['orbital_energy'] for d in data]
            z0,z1=[complex(*d[key]) for d in data]
            if not np.isfinite([before,after,z0,z1]).all():raise ValueError('Nonfinite flux/amplitude')
            relative=(after-before)/before if before else None
            if before:
                np.testing.assert_allclose(relative,abs(z1/z0)**2-1,atol=2e-14,rtol=1e-9)
            if boundary=='infinity' and not data[0]['propagating'] and (before!=0 or after!=0):
                raise ValueError('Bound mode has nonzero infinity flux')
            flux[boundary]=dict(default=before,dense=after,absolute_change=after-before,relative_change=relative)
        norm=np.dot(weights,abs(source[0])**2)
        rows.append(dict(ell=mode[0],m=mode[1],propagating=data[0]['propagating'],samples=len(weights),
            inputs=[dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths],
            source_weighted_relative_l2=float(np.sqrt(np.dot(weights,abs(source[1]-source[0])**2)/norm)) if norm else None,
            flux=flux))
    return dict(status='completed_backend_comparison_not_full_convergence',
        inputs=[dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in (default_path,dense_path)],rows=rows,
        limitation='Same finite grid and truncation. Source differences are discrete norms, not continuum error bounds.')


def main():
    p=argparse.ArgumentParser();p.add_argument('default',type=Path);p.add_argument('dense',type=Path);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();result=compare(args.default,args.dense)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    for row in result['rows']:print(row['ell'],row['m'],row['flux'])


if __name__=='__main__':main()
