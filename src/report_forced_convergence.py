"""Compare completed actual-source runs without asserting global convergence."""
import json
from pathlib import Path


def main():
    directory=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    rows=[]
    previous=None
    physical=None
    for order in (2,4,8,16,32):
        path=directory/f'forced_mode_nr{order}_nt6_L4.json'
        if not path.exists():
            continue
        data=json.loads(path.read_text())
        if data['status']!='truncated_single_mode_not_converged':
            continue
        parameters={k:v for k,v in data['parameters'].items() if k!='radial_order'}
        if physical is not None and parameters!=physical:
            raise ValueError('Cannot compare runs with different physics or cutoffs')
        physical=parameters
        row=dict(radial_order=order,sample_count=len(data['samples']),
                 infinity_orbital_flux=data['flux']['infinity']['orbital_energy'],
                 horizon_orbital_flux=data['flux']['horizon']['orbital_energy'])
        if previous:
            for key in ('infinity_orbital_flux','horizon_orbital_flux'):
                row[key+'_relative_change']=abs(row[key]-previous[key])/max(abs(row[key]),1e-300)
        if 'outer_source_cutoff_sequence' in data:
            before,after=data['outer_source_cutoff_sequence'][-2:]
            row['last_outer_cutoffs']=[before['outer_source_cutoff'],after['outer_source_cutoff']]
            for key in ('z_inf','z_h'):
                a,b=complex(*before[key]),complex(*after[key])
                row[key+'_outer_tail_relative_change']=abs(a-b)/max(abs(b),1e-300)
        rows.append(row)
        previous=row
    log_rows=[]
    for order in (8,16,32,64):
        path=directory/f'forced_mode_nr8_nt6_L4_log_h{order}.json'
        if not path.exists():
            continue
        data=json.loads(path.read_text())
        if data['status']!='truncated_single_mode_not_converged':
            continue
        row=dict(horizon_order=order,horizon_orbital_flux=data['flux']['horizon']['orbital_energy'],
                 infinity_orbital_flux=data['flux']['infinity']['orbital_energy'],
                 reused_source_samples=sum(sample.get('reused',False) for sample in data['samples']))
        if log_rows:
            row['horizon_flux_relative_change']=abs(row['horizon_orbital_flux']-log_rows[-1]['horizon_orbital_flux'])/max(abs(row['horizon_orbital_flux']),1e-300)
        log_rows.append(row)
    cutoff_rows=[]
    for offset in (.05,.005,.0005):
        name=('forced_mode_nr8_nt6_L4_log_h32.json' if offset==.05 else
              f'forced_mode_nr8_nt6_L4_inner{offset:g}_outer320_log_h32.json')
        path=directory/name
        if not path.exists():
            continue
        data=json.loads(path.read_text())
        if data['status']!='truncated_single_mode_not_converged':
            continue
        row=dict(source_inner_offset=offset)
        for boundary in ('infinity','horizon'):
            key=boundary+'_orbital_flux'
            row[key]=data['flux'][boundary]['orbital_energy']
            if cutoff_rows:
                row[key+'_relative_change']=abs(row[key]-cutoff_rows[-1][key])/max(abs(row[key]),1e-300)
        cutoff_rows.append(row)
    refinement={}
    for label,settings in (
        ('angular',[(n,4) for n in (6,10,14)]),
        ('metric',[(10,l) for l in (4,6,8,12,18)])):
        sequence=[]
        for n,l in settings:
            path=directory/f'forced_mode_nr8_nt{n}_L{l}_log_h16.json'
            if not path.exists():
                continue
            data=json.loads(path.read_text())
            if data['status']!='truncated_single_mode_not_converged':
                continue
            row=dict(angular_order=n,metric_ellmax=l)
            for boundary in ('infinity','horizon'):
                key=boundary+'_orbital_flux'
                row[key]=data['flux'][boundary]['orbital_energy']
                if sequence:
                    row[key+'_relative_change']=abs(row[key]-sequence[-1][key])/max(abs(row[key]),1e-300)
            sequence.append(row)
        refinement[label]=sequence
    result=dict(status='diagnostic_not_full_paper_convergence',runs=rows,
                horizon_log_runs=log_rows,source_inner_cutoff_runs=cutoff_rows,
                angular_refinement=refinement['angular'],metric_refinement=refinement['metric'])
    markers=directory/'paper_v1_flux_markers.json'
    if markers.exists():
        reference=json.loads(markers.read_text())
        plotted=next(row['plotted_flux'] for row in reference['markers'] if row['ell']==3 and row['m']==3)
        comparisons=[]
        for row in refinement['metric']:
            value=row['infinity_orbital_flux']
            comparisons.append(dict(angular_order=row['angular_order'],metric_ellmax=row['metric_ellmax'],
                flux_per_q2_eta=value,flux_per_q2_epsilon2_from_text=value/.3**6,
                ratio_unit_mass_to_plotted=value/plotted,
                ratio_text_epsilon_to_plotted=value/(.3**6*plotted)))
        result['paper_figure6_normalization_audit']=dict(status='unresolved_not_a_reproduction_error_estimate',
            plotted_l3_m3=plotted,epsilon2_over_eta=.3**6,
            reference='paper_v1_flux_markers.json; digitized graphic, not author numeric data',
            comparisons=comparisons)
    (directory/'forced_mode_convergence.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
