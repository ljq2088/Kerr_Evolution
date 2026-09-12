"""Relate saved actual scalar (0,0) responses to the independently found pole."""
import hashlib
import json
from pathlib import Path
import numpy as np


def main():
    folder=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    pole_path=folder/'kerr_resonance_poles.json'
    spectrum=json.loads(pole_path.read_text())
    case=next(c for c in spectrum['cases'] if (c['ell'],c['m'])==(0,0))
    pole=complex(*case['refined']['omega']);width=-pole.imag
    if width<=0:raise ValueError('Expected a damped |100> pole')
    paths=[folder/'flux_coverage_L18_nt18_i6_h5_f5_rp10.json',
           folder/'flux_coverage_L18_nt18_i6_h5_f5_m0go4000_m2go4000.json',
           folder/'flux_coverage_L18_nt18_i6_h5_f5_rp30.json']
    for radius in (41.6,41.8):
        paths.append(next(folder.glob(f'flux_coverage*_m1go4000_*rp{radius}.json')))
    inputs={pole_path.name:hashlib.sha256(pole_path.read_bytes()).hexdigest()};rows=[]
    for coverage_path in paths:
        coverage=json.loads(coverage_path.read_text())
        candidates=[r for r in coverage['horizon']['modes'] if (r['ell'],r['m'])==(0,0)]
        if len(candidates)!=1 or 'flux' not in candidates[0]:raise ValueError('Missing (0,0) response')
        path=folder/candidates[0]['file'];d=json.loads(path.read_text());p=d['parameters'];metric=p['metric']
        if d['status']!='truncated_single_mode_not_converged' or p.get('background') is not None:
            raise ValueError('Require a completed stationary Kerr response')
        if p['alpha']!=spectrum['alpha'] or metric['a']!=spectrum['a'] or metric['m_g']!=-1:
            raise ValueError('Pole and forcing backgrounds differ')
        radius=metric['orbital_radius'];op=1/(radius**1.5+metric['a'])
        omega=p['omega'];zh=complex(*d['z_h']);horizon=1+np.sqrt(1-metric['a']**2)
        if abs(omega+op-spectrum['cloud_frequency'])>1e-13:raise ValueError('Wrong drive frequency')
        flux=d['flux']['horizon']['orbital_energy']
        np.testing.assert_allclose(flux,-4*horizon*op*omega*abs(zh)**2,rtol=1e-12,atol=1e-30)
        inputs[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        detuning=omega-pole.real
        numerator=(omega-pole)*zh
        rows.append(dict(rp=radius,drive_frequency=omega,detuning=detuning,
            detuning_over_width=detuning/width,horizon_orbital_flux=flux,
            z_h=[zh.real,zh.imag],z_h_magnitude=abs(zh),
            pole_removed_response_magnitude=abs(numerator),
            green_outer_radius=p['green_outer_radius'],infinity_method=p.get('infinity_method','series')))
    result=dict(status='actual_forced_samples_with_pole_not_resonance_scan',
        alpha=spectrum['alpha'],a=spectrum['a'],pole=[pole.real,pole.imag],
        real_frequency_crossing_radius=case['resonance_radius_from_real_frequency'],
        rows=rows,inputs=inputs,amplitude_fit=False,
        limitation='The pole-removed quantity includes varying source coupling and regular Green contributions; it is not a measured pole residue. Sparse finite-resolution samples do not locate a flux maximum or prove a floating orbit.')
    (folder/'forced_100_pole_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')
    for row in rows:print(row['rp'],row['detuning_over_width'],row['horizon_orbital_flux'],row['pole_removed_response_magnitude'])


if __name__=='__main__':main()
