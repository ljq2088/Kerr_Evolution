"""Compare physical static responses at identical source nodes and panel radii."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('before', type=Path)
    parser.add_argument('after', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    batches = [json.loads(p.read_text()) for p in (args.before, args.after)]
    for batch in batches:
        if batch['status'] != 'batch_completed_finite_resolution_not_converged':
            raise ValueError('Require completed batches')
    if batches[0]['metric_precomputation']['source_hash'] != batches[1]['metric_precomputation']['source_hash']:
        raise ValueError('Different metric source versions')
    maps = [{(r['scalar_ell'], r['scalar_m']): r['file'] for r in b['channels']} for b in batches]
    if maps[0].keys() != maps[1].keys():
        raise ValueError('Different channel sets')
    rows = []
    for key in sorted(maps[0]):
        paths = [p.parent / m[key] for p, m in zip((args.before, args.after), maps)]
        records = [json.loads(p.read_text()) for p in paths]
        parameters = [dict(d['parameters']) for d in records]
        boundaries = [p.pop('green_outer_radius') for p in parameters]
        if parameters[0] != parameters[1] or key[1] != 1:
            raise ValueError('Only the static scalar Green outer boundary may differ')
        samples = [[(s['r'], s['weight'], s['source']) for s in d['samples']] for d in records]
        if samples[0] != samples[1]:
            raise ValueError('Source nodes, weights or values differ')
        for d in records:
            if d['status'] != 'truncated_single_mode_not_converged' or any(
                v != 0 for boundary in d['flux'].values() for v in boundary.values()):
                raise ValueError('Unexpected static response status or nonzero threshold flux')
        fields = [{r['r']: complex(*r['field']) for r in d['radial_response_at_panel_boundaries']} for d in records]
        radii = sorted(set(fields[0]) & set(fields[1]))
        old, new = [np.array([f[r] for r in radii]) for f in fields]
        r0 = parameters[0]['metric']['orbital_radius']
        point_old, point_new = [f[r0] for f in fields]
        zh_old, zh_new = [complex(*d['z_h']) for d in records]
        rows.append(dict(ell=key[0], m=key[1], outer_boundaries=boundaries,
            inputs={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
            identical_source_samples=len(samples[0]), common_field_radii=radii,
            particle_field_relative_change=abs((point_new-point_old)/point_old),
            horizon_amplitude_relative_change=abs((zh_new-zh_old)/zh_old),
            maximum_absolute_field_change=float(max(abs(new-old))),
            max_field_change_over_max_field=float(max(abs(new-old))/max(abs(old)))))
    result = dict(status='static_same_source_boundary_check_not_full_convergence', rows=rows,
        limitations=['Physical fields compared; arbitrary bound Up normalization is not compared',
            'Only sampled common radii checked; not a uniform continuum error bound',
            'Source discretization, static matching and metric cutoff remain fixed'])
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(rows), flush=True)


if __name__ == '__main__':
    main()
