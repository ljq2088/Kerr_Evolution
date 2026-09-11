"""Single-channel Kerr/Schwarzschild comparison, explicitly not total flux."""
import json
from pathlib import Path


def main():
    folder=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    files=dict(kerr='forced_mode_nr8_nt10_L4_log_h16.json',
               schwarzschild='forced_mode_nr8_nt10_L4_schwarzschild_frozen_log_h16.json')
    runs={}
    for label,name in files.items():
        data=json.loads((folder/name).read_text())
        if data['status']!='truncated_single_mode_not_converged':
            raise ValueError('Comparison requires completed integrals')
        p=data['parameters']
        for key,value in dict(alpha=.3,radial_order=8,angular_order=10,scalar_ell=3,scalar_m=3).items():
            if p[key]!=value:raise ValueError(f'Unexpected comparison parameter {key}')
        if p['metric']['ellmax']!=4 or p['metric']['orbital_radius']!=20.:
            raise ValueError('Mismatched orbit/metric truncation')
        runs[label]=dict(file=name,parameters=p,flux=data['flux'])
    s=runs['schwarzschild']['flux']['infinity']['orbital_energy']
    k=runs['kerr']['flux']['infinity']['orbital_energy']
    result=dict(status='single_mode_diagnostic_not_figure2_or_3',
        runs=runs,infinity_relative_difference_to_schwarzschild=(k-s)/s,
        limitations=['Only scalar (ell,m)=(3,3), not a full flux sum',
                     'Schwarzschild temporal decay is frozen; complex radial cloud retained',
                     'Unit cloud Killing energy uses the documented slicing conventions',
                     'Source and boundary convergence incomplete; paper normalization unresolved'])
    (folder/'background_single_mode_comparison.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    main()
