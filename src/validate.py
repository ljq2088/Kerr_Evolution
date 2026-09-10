"""Reproducible radial/angular/time convergence and Kerr waveform examples."""
from dataclasses import replace
from pathlib import Path
import json
import numpy as np
from kerr_scalar import Config, KerrScalar, save_run


def waves(config):
    _, f, _ = KerrScalar(config).evolve()
    return f[:, 0, [0, -1], :]


def main():
    output = Path('outputs/validation')
    output.mkdir(parents=True, exist_ok=True)
    report = {}
    for m in (0, 2):
        c = Config(m=m, lmax=6, tmax=80, samples=201, rtol=2e-11, atol=2e-13)
        radial = [waves(replace(c, n=n)) for n in (32, 48, 64)]
        differences = [float(np.max(np.abs(radial[i]-radial[i+1]))) for i in (0, 1)]
        angular4 = waves(replace(c, n=64, lmax=4))
        angular8 = waves(replace(c, n=64, lmax=8))
        e46 = float(np.max(np.abs(angular4-radial[-1][:, :, :5-m])))
        e68 = float(np.max(np.abs(radial[-1]-angular8[:, :, :7-m])))
        temporal = waves(replace(c, n=64, rtol=1e-9, atol=1e-11))
        et = float(np.max(np.abs(temporal-radial[-1])))
        report[f'm={m}'] = dict(radial_n=[32,48,64], radial_max_abs_differences=differences,
                               angular_lmax=[4,6,8], angular_common_modes_max_abs_differences=[e46,e68],
                               temporal_max_abs_difference=et,
                               peak_boundary_amplitude=float(np.max(np.abs(radial[-1]))))
        assert differences[1] < differences[0]/10, report
        assert differences[1] < 2e-5, report
        assert e68 < e46/5, report
        assert et < 1e-7, report
    (output/'convergence.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)
    save_run(Config(n=64, lmax=6, m=0), 'outputs/kerr_m0')
    save_run(Config(n=64, lmax=6, m=2), 'outputs/kerr_m2')


if __name__ == '__main__':
    main()
