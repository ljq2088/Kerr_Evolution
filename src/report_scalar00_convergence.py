"""Actual-source scalar (0,0) channel, separate one-parameter convergence tests."""
import json
from pathlib import Path
import numpy as np


def main():
    folder=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    files=dict(baseline='forced_mode_nr4_nt10_L4_mg-1_sl0_log_h16.json',
        radial='forced_mode_nr8_nt10_L4_mg-1_sl0_log_h16.json',
        metric='forced_mode_nr4_nt10_L6_mg-1_sl0_log_h16.json',
        inner005='forced_mode_nr4_nt10_L4_mg-1_sl0_inner0.005_outer320_log_h16.json',
        inner0005='forced_mode_nr4_nt10_L4_mg-1_sl0_inner0.0005_outer320_log_h16.json',
        inner0005_h32='forced_mode_nr4_nt10_L4_mg-1_sl0_inner0.0005_outer320_log_h32.json')
    cases={};pending=[]
    for label,name in files.items():
        path=folder/name
        if not path.exists():pending.append(name);continue
        data=json.loads(path.read_text())
        if 'flux' not in data:pending.append(name);continue
        p=data['parameters'];m=p['metric']
        if (p['alpha'],p['scalar_ell'],p['scalar_m'],m['m_g'],m['orbital_radius'])!=(.3,0,0,-1,20.):
            raise ValueError('Unexpected physical channel')
        a=m['a'];rp=1+np.sqrt(1-a*a)
        cases[label]=dict(file=name,radial_order=p['radial_order'],angular_order=p['angular_order'],
            horizon_order=p['horizon_quadrature_order'],
            metric_ellmax=m['ellmax'],source_inner_offset=p['source_panels'][0]-rp,
            infinity_flux=data['flux']['infinity']['orbital_energy'],
            horizon_flux=data['flux']['horizon']['orbital_energy'],
            horizon_wave_energy=data['flux']['horizon']['wave_energy'],
            horizon_charge=data['flux']['horizon']['charge'],omega=p['omega'])
    baseline=cases['baseline']
    for name in ('radial','metric','inner005','inner0005','inner0005_h32'):
        if name in cases:
            cases[name]['relative_change_from_baseline']=abs(cases[name]['horizon_flux']/baseline['horizon_flux']-1)
    if 'inner0005_h32' in cases:
        cases['inner0005_h32']['relative_change_from_same_cutoff_h16']=abs(cases['inner0005_h32']['horizon_flux']/cases['inner0005']['horizon_flux']-1)
    a=json.loads((folder/files['baseline']).read_text())['parameters']['metric']['a']
    Omega=1/(20**1.5+a);omega_c=baseline['omega']+Omega
    loss=-omega_c*baseline['horizon_charge']
    balance=baseline['horizon_wave_energy']+loss-baseline['horizon_flux']
    result=dict(status='single_channel_convergence_incomplete',flux_scaling='per q^2 eta, eta=M_cloud/M; not epsilon-normalized',
        physical_mode='scalar ell=m=0 is sourced by nonstatic metric m_g=-1',cases=cases,pending=pending,
        orbit_frequency=Omega,cloud_frequency=omega_c,
        baseline_cloud_energy_change=loss,baseline_energy_accounting_residual=balance,
        interpretation='Positive scalar energy enters the horizon; cloud depletion yields a negative orbital-effective flux. No full-flux or floating-orbit claim.',
        missing=['near-horizon quadrature and cutoff convergence','angular projection check','full mode sums','paper normalization reconciliation'])
    (folder/'scalar00_convergence.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    main()
