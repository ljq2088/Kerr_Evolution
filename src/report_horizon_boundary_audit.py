"""Refine only the radial Green boundary/tolerance for a saved source."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from environment_radial import RadialGreen


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--outer', type=float, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    if data.get('status') != 'truncated_single_mode_not_converged':
        raise ValueError('Require a completed finite response')
    if data.get('flux_validity') == 'historical_unreliable_boundary_result':
        raise ValueError('Invalid historical boundary result')
    p = data['parameters']; metric = p['metric']
    original_outer = p['green_outer_radius']
    if args.outer <= original_outer:
        raise ValueError('Refined boundary must lie outside the original boundary')
    r = np.array([s['r'] for s in data['samples']])
    w = np.array([s['weight'] for s in data['samples']])
    source = np.array([complex(*s['source']) for s in data['samples']])
    recorded = complex(*data['z_h'])
    if not recorded:
        raise ValueError('Cannot form relative changes for a zero response')
    runs = []
    for outer, tolerance in ((original_outer, 1e-11), (args.outer, 1e-11),
                             (args.outer, 1e-13)):
        green = RadialGreen(metric['a'], p['alpha'], p['omega'],
            p['scalar_ell'], p['scalar_m'], rmax=outer,
            offset=p['green_horizon_offset'], rtol=tolerance,
            infinity_method=p.get('infinity_method', 'series'))
        terms = w * (green.upsol.sol(r)[0] / green.w0) * source
        amplitude = terms.sum()
        if not np.isfinite(amplitude):
            raise ValueError('Nonfinite response')
        change = abs(amplitude / recorded)**2 - 1
        row = dict(outer=outer, rtol=tolerance,
            z_h=[float(amplitude.real), float(amplitude.imag)],
            relative_complex_amplitude_change=float(abs(amplitude/recorded-1)),
            relative_horizon_flux_change=float(change),
            horizon_orbital_energy=data['flux']['horizon']['orbital_energy']*(1+change),
            discrete_condition_number=float(sum(abs(terms))/abs(amplitude)),
            wronskian_relative_spread=float(max(abs(green.wronskian(r)/green.w0-1))),
            last_series_term=green.series_last_term_relative)
        if not runs and row['relative_complex_amplitude_change'] > 1e-9:
            raise ValueError('Baseline does not reproduce saved amplitude')
        runs.append(row); print(json.dumps(row), flush=True)
    result = dict(status='same_source_horizon_boundary_audit_not_full_convergence',
        input=args.input.name, sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
        parameters=p, runs=runs,
        limitations=['Source nodes and weights remain fixed',
            'No source cutoff, source quadrature, metric or angular convergence test',
            'Discrete conditioning is not a continuum error bound'])
    args.output.write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
