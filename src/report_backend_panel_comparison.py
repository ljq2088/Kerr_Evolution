"""Locate backend-induced amplitude changes in saved radial source panels.

Uses the existing cumulative integrals, without recomputing the metric or Green
function. This is discrete attribution, not an estimate of physical error.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from report_backend_batch_comparison import compare


def encode(z):
    return [float(z.real), float(z.imag)]


def panel_changes(before, after):
    edges = np.asarray(before['parameters']['source_panels'], dtype=float)
    if not np.array_equal(edges, after['parameters']['source_panels']):
        raise ValueError('Different source panels')
    if not np.isfinite(edges).all() or np.any(np.diff(edges) <= 0):
        raise ValueError('Invalid source panels')
    seq = [d['outer_source_cutoff_sequence'] for d in (before, after)]
    for rows in seq:
        if [r['outer_source_cutoff'] for r in rows] != edges[1:].tolist():
            raise ValueError('Incomplete cumulative panel integrals')
    result = {}
    for key in ('z_h', 'up_coefficient'):
        cumulative = [np.array([0j] + [complex(*r[key]) for r in rows]) for rows in seq]
        if not np.isfinite(cumulative).all():
            raise ValueError('Nonfinite cumulative amplitude')
        for d, values in zip((before, after), cumulative):
            np.testing.assert_allclose(values[-1], complex(*d[key]), rtol=2e-13, atol=0)
        panels = [np.diff(values) for values in cumulative]
        change = panels[1] - panels[0]
        total = cumulative[1][-1] - cumulative[0][-1]
        scale = max(abs(cumulative[0][-1]), abs(cumulative[1][-1]))
        np.testing.assert_allclose(change.sum(), total, rtol=1e-10, atol=2e-14*scale)
        norm = float(np.sum(abs(change)))
        largest = int(np.argmax(abs(change)))
        result[key] = dict(
            total_change=encode(total), sum_absolute_panel_changes=norm,
            total_change_over_sum_absolute=float(abs(total)/norm) if norm else None,
            largest_change_interval=edges[largest:largest+2].tolist(),
            panels=[dict(inner=float(lo), outer=float(hi), default=encode(z0),
                         dense=encode(z1), change=encode(dz))
                    for lo,hi,z0,z1,dz in zip(edges[:-1],edges[1:],*panels,change)])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('default', type=Path)
    parser.add_argument('dense', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = compare(args.default, args.dense)
    for row in result['rows']:
        paths = [p.parent / item['file'] for p,item in
                 zip((args.default,args.dense),row['inputs'])]
        row['radial_panel_amplitudes'] = panel_changes(*[json.loads(p.read_text()) for p in paths])
    result['limitation'] += ' Panel changes attribute the saved finite sum only; cancellation may suppress the total.'
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(args.output)


if __name__ == '__main__':
    main()
