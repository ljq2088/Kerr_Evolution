"""Export the legacy scalar solution in the upstream Ripley coordinates."""
import argparse
import json
from pathlib import Path
import numpy as np


def convert(source, destination, length_scale=None):
    with np.load(source, allow_pickle=False) as data:
        config = json.loads(str(data['config']))
        mass = config['mass']
        length = mass if length_scale is None else float(length_scale)
        if not np.isfinite(length) or length <= 0:
            raise ValueError('length_scale must be finite and positive')
        rp = mass * (1 + np.sqrt(1 - config['spin']**2))
        offset = 4 * mass * np.log(rp / mass)
        metadata = dict(config=config, spin_weight=0, L=length, r_plus=rp,
                        tau_minus_T=offset, field='Phi_scalar = R * psi0',
                        azimuth='ingoing regular Phi; exp(i*m*Phi)',
                        discretization='Chebyshev / spherical Galerkin / DOP853')
        destination = Path(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(destination, T=data['tau']-offset,
                            T_elapsed=data['tau']-data['tau'][0],
                            R=length**2/rp*data['sigma'], ell=data['ell'],
                            psi0=data['u']/length**2,
                            psi0_T=data['p']/length**2,
                            metadata=json.dumps(metadata))
    return metadata


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source')
    parser.add_argument('destination')
    parser.add_argument('--length-scale', type=float)
    args = parser.parse_args()
    print(json.dumps(convert(args.source, args.destination, args.length_scale), indent=2))
