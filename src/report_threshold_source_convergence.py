"""Track actual near-threshold source refinement separately from boundary tests."""
import json
from pathlib import Path


def main():
    folder=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    runs=[]
    for path in folder.glob('forced_mode_*_alpha0.3_rp41.6_mg1_sl2_gh0.0001_go32000_coulomb_log_h*.json'):
        data=json.loads(path.read_text())
        if data['status']!='truncated_single_mode_not_converged':
            continue
        p=data['parameters']
        runs.append(dict(file=path.name,radial_order=p['radial_order'],angular_order=p['angular_order'],
            metric_ellmax=p['metric']['ellmax'],horizon_order=p['horizon_quadrature_order'],
            infinity_flux=data['flux']['infinity']['orbital_energy'],
            horizon_flux=data['flux']['horizon']['orbital_energy']))
    runs.sort(key=lambda row:(row['metric_ellmax'],row['angular_order'],row['radial_order']))
    for previous,row in zip(runs[:-1],runs[1:]):
        row['changed_resolution_parameters']=[key for key in ('radial_order','angular_order','metric_ellmax','horizon_order')
                                               if row[key]!=previous[key]]
        for key in ('infinity_flux','horizon_flux'):
            row[key+'_relative_change']=abs(row[key]-previous[key])/max(abs(row[key]),1e-300)
    (folder/'threshold_source_convergence.json').write_text(json.dumps(dict(
        status='not_converged_paper_result',orbit=41.6,alpha=.3,scalar_mode=[2,2],
        flux_scaling='per q² eta; paper epsilon normalization unresolved',runs=runs),indent=2)+'\n')


if __name__=='__main__':
    main()
