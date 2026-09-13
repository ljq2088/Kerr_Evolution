"""Coherent dipole/higher-degree attribution, not a sum of separate fluxes."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import numpy as np


def compare(full_path,dipole_path):
    full,dipole=[json.loads(p.read_text()) for p in (full_path,dipole_path)]
    for d in (full,dipole):
        if d.get('status')!='truncated_single_mode_not_converged':raise ValueError('Both responses must be complete')
    params=[copy.deepcopy(d['parameters']) for d in (full,dipole)]
    limits=[p['metric'].pop('ellmax') for p in params]
    if limits!=[18,1] or params[0]!=params[1]:raise ValueError('Only metric ellmax18 versus1 may differ')
    if [[(r['r'],r['weight']) for r in d['samples']] for d in (full,dipole)][0]!=[(r['r'],r['weight']) for r in dipole['samples']]:
        raise ValueError('Source quadrature differs')
    rows={}
    for boundary,key in [('horizon','z_h'),('infinity','z_inf')]:
        z,d=[complex(*x[key]) for x in (full,dipole)];other=z-d
        flux=full['flux'][boundary]['orbital_energy']
        if not np.isfinite([z,d,flux]).all():raise ValueError('Nonfinite amplitude')
        if not z:
            if flux!=0 or dipole['flux'][boundary]['orbital_energy']!=0:raise ValueError('Inconsistent zero amplitude')
            rows[boundary]=dict(full_flux=0.,dipole_flux=0.,higher_degree_flux=0.,interference_flux=0.)
            continue
        coefficient=flux/abs(z)**2
        alone=coefficient*abs(d)**2;higher=coefficient*abs(other)**2;cross=2*coefficient*(d*other.conjugate()).real
        np.testing.assert_allclose(alone,dipole['flux'][boundary]['orbital_energy'],rtol=2e-12,atol=1e-30)
        np.testing.assert_allclose(alone+higher+cross,flux,rtol=2e-12,atol=1e-30)
        rows[boundary]=dict(full_flux=flux,dipole_flux=alone,higher_degree_flux=higher,interference_flux=cross,
            dipole_amplitude_relative_to_full=[(d/z).real,(d/z).imag],
            higher_degree_amplitude_relative_to_full=[(other/z).real,(other/z).imag])
    return dict(status='coherent_metric_degree_attribution_not_error_estimate',
        inputs=[dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in (full_path,dipole_path)],
        scalar_mode=[params[0]['scalar_ell'],params[0]['scalar_m']],rows=rows,
        limitation='Metric separated degree, not spherical degree. Nonzero contributions alone do not identify a faulty sector. Full flux includes interference.')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('full',type=Path);p.add_argument('dipole',type=Path);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();result=compare(args.full,args.dipole);args.output.write_text(json.dumps(result,indent=2)+'\n');print(result['rows'])


if __name__=='__main__':main()
