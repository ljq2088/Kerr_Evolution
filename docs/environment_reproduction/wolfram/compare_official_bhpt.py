"""Compare untouched official Mathematica BHPT samples against installed pybhpt AUTO.
Run from repository root with .venv/bin/python docs/environment_reproduction/wolfram/compare_official_bhpt.py.
Official scalar/complex results are read only; no phase or amplitude fitting is needed for the raw comparison.
"""
from pathlib import Path
import json
import numpy as np
from pybhpt.radial import RadialTeukolsky
from pybhpt.swsh import SpinWeightedSpheroidalHarmonic

ROOT = Path(__file__).resolve().parents[3]
SAMPLES = ROOT / 'outputs/paper_metric_reference/wolfram_dependencies/official_results'
def complex_value(x):
    return complex(x['re'], x['im']) if isinstance(x, dict) else complex(x)
def complex_json(x):
    return {'re': float(x.real), 'im': float(x.imag)}

def main():
    results = []
    for path in sorted(SAMPLES.glob('*_s*_m*.json')):
        record = json.loads(path.read_text()); assert record['ok']
        v = record['data']
        s, ell, m = map(int, (v['s'], v['l'], v['m']))
        a, omega = v['a'], v['omega']
        radii = np.array([row['r'] for row in v['radial']['In']['samples']])
        rt = RadialTeukolsky(s, ell, m, a, omega, radii); rt.solve(method='AUTO')
        sw = SpinWeightedSpheroidalHarmonic(s, ell, m, a * omega)
        theta = np.array([row['theta'] for row in v['angular']])
        angular = np.array([complex_value(row['S']) for row in v['angular']])
        angular_derivative = np.array([complex_value(row['dS']) for row in v['angular']])
        angular_py = sw.Sslm(theta); angular_derivative_py = sw.Sslm_derivative(theta)
        assert np.all(np.isfinite([angular, angular_derivative, angular_py, angular_derivative_py]))
        boundaries, original, native = {}, {}, {}
        for bc in ('In', 'Up'):
            rows = v['radial'][bc]['samples']
            field = np.array([complex_value(row['R']) for row in rows])
            derivative = np.array([complex_value(row['dR']) for row in rows])
            field_py, derivative_py = rt.radialsolutions(bc), rt.radialderivatives(bc)
            assert np.all(np.isfinite([field, derivative, field_py, derivative_py]))
            ratio = field_py[1] / field[1]
            boundaries[bc] = {
                'R_rel_max': float(np.max(np.abs((field_py - field) / field))),
                'dR_rel_max': float(np.max(np.abs((derivative_py - derivative) / derivative))),
                'rescaled_R_rel_max': float(np.max(np.abs((field_py / ratio - field) / field))),
                'rescaled_dR_rel_max': float(np.max(np.abs((derivative_py / ratio - derivative) / derivative))),
                'normalization_ratio_at_middle': complex_json(ratio),
                'official_transmission': v['radial'][bc]['amplitudes']['Transmission'],
            }
            original[bc] = field, derivative; native[bc] = field_py, derivative_py
        delta = radii**2 - 2 * radii + a*a
        wronskian = (original['In'][0] * original['Up'][1] - original['Up'][0] * original['In'][1]) * delta**(s+1)
        wronskian_py = (native['In'][0] * native['Up'][1] - native['Up'][0] * native['In'][1]) * delta**(s+1)
        kernel = original['In'][0][0] * original['Up'][0][-1] / wronskian[1]
        kernel_py = native['In'][0][0] * native['Up'][0][-1] / wronskian_py[1]
        angular_anchor = int(np.argmax(np.abs(angular)))  # avoid exact equatorial nodes
        result = {
            'case': v['case'], 'angular_eigen_abs': float(abs(sw.eigenvalue - v['eigenvalue'])),
            'angular_S_abs_max': float(np.max(np.abs(angular_py - angular))),
            'angular_dS_abs_max': float(np.max(np.abs(angular_derivative_py - angular_derivative))),
            'angular_ratio_at_max_sample': complex_json(angular_py[angular_anchor] / angular[angular_anchor]),
            'radial': boundaries,
            'official_W_drift': float(np.max(np.abs(wronskian / wronskian[1] - 1))),
            'native_W_drift': float(np.max(np.abs(wronskian_py / wronskian_py[1] - 1))),
            'physical_G_rel': float(abs(kernel_py / kernel - 1)),
        }
        results.append(result)
    assert len(results) == 11
    summary = {
        'pybhpt': '1.0.0', 'method': 'AUTO', 'case_count': len(results),
        'kernel_definition': 'R_In(r_min) R_Up(r_max) / [Delta^(s+1) (R_In R_Up_prime - R_Up R_In_prime)] evaluated at middle radius',
        'maximum_R_relative': max(x['radial'][bc]['R_rel_max'] for x in results for bc in ('In','Up')),
        'maximum_dR_relative': max(x['radial'][bc]['dR_rel_max'] for x in results for bc in ('In','Up')),
        'maximum_angular_eigen_absolute': max(x['angular_eigen_abs'] for x in results),
        'maximum_angular_S_absolute': max(x['angular_S_abs_max'] for x in results),
        'maximum_angular_dS_absolute': max(x['angular_dS_abs_max'] for x in results),
        'maximum_official_W_drift': max(x['official_W_drift'] for x in results),
        'maximum_physical_G_relative': max(x['physical_G_rel'] for x in results),
        'cases': results,
    }
    encoded = json.dumps(summary, indent=2, allow_nan=False) + '\n'
    (SAMPLES / 'pybhpt_AUTO_comparison.json').write_text(encoded)
    print(json.dumps({k:v for k,v in summary.items() if k != 'cases'}, indent=2))

if __name__ == '__main__':
    main()
