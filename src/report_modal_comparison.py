"""Plot complete ell<=6 finite-resolution inventory against v1 markers.

This is a diagnostic; it does not include Fig.6's ell=7..12 or claim that
the metric, source, boundary, or normalization has converged.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--metric-ellmax',type=int,default=6)
    parser.add_argument('--coverage',type=Path,help='Explicit completed ell<=6 / horizon ell<=5 inventory')
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    coverage_path=args.coverage or root/f'flux_coverage_L{args.metric_ellmax}.json'
    coverage=json.loads(coverage_path.read_text())
    cutoff=coverage['parameters']['metric_ellmax']
    if (coverage['infinity']['required_count']!=9 or coverage['horizon']['required_count']!=18
        or any(row['ell']>6 for row in coverage['infinity']['modes'])
        or any(row['ell']>5 for row in coverage['horizon']['modes'])):
        raise ValueError('This comparison requires infinity ell<=6 and horizon ell<=5')
    if any(coverage[b]['finite_resolution_total'] is None for b in ('infinity','horizon')):
        raise ValueError('Complete the common-resolution flux inventory first')
    markers=json.loads((root/'paper_v1_flux_markers.json').read_text())
    figure2=json.loads((root/'paper_normalization_audit.json').read_text())
    rows=[]
    for row in coverage['infinity']['modes']:
        reference=next(r for r in markers['markers'] if (r['ell'],r['m'])==(row['ell'],row['m']))
        actual=row['flux']['infinity']['orbital_energy'];target=reference['plotted_flux']
        rows.append(dict(ell=row['ell'],m=row['m'],file=row['file'],computed=actual,
                         digitized=target,relative_difference=actual/target-1))
    report=dict(status='finite_resolution_comparison_not_reproduction',parameters=coverage['parameters'],
        coverage_sha256=hashlib.sha256(coverage_path.read_bytes()).hexdigest(),
        modes=rows,finite_infinity_total=coverage['infinity']['finite_resolution_total'],
        finite_horizon_total=coverage['horizon']['finite_resolution_total'],
        sum_figure6_ell_le6=sum(row['digitized'] for row in rows),
        figure2_horizon_magnitude_inferred_unit_cloud_mass=figure2['figure2']['horizon_magnitude']['plotted_flux']*.3**6,
        limitations=['Raw finite source cutoffs, no horizon-tail correction',
                     ('Metric ell cutoff below paper target 18' if cutoff<18 else 'Reconstruction index cutoff 18 reached; equivalence to paper spherical-data truncation unverified'),
                     'Scalar ell=7..12 of paper Fig.6 not present',
                     'Figure 2/6 normalization relationship inferred from cross-figure closure'])
    stem=root/f'rp20_modal_comparison_L{cutoff}'
    stem.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    x=np.arange(len(rows))
    fig,(ax,err)=plt.subplots(2,1,figsize=(9,6),sharex=True,layout='constrained',height_ratios=[2,1])
    ax.semilogy(x,[r['digitized'] for r in rows],'o',color='.4',label='Paper v1 vector markers')
    ax.semilogy(x,[r['computed'] for r in rows],'x',color='#1764ab',ms=8,label=f'Computed metric L={cutoff}')
    ax.set(ylabel=r'Infinity orbital flux per $q^2(M_c/M)$',title='Kerr, alpha=0.3, r_p=20M — finite-resolution diagnostic')
    ax.legend();ax.grid(alpha=.2)
    err.bar(x,[100*r['relative_difference'] for r in rows],color='#1764ab')
    err.axhline(0,color='.3',lw=.8);err.set(ylabel='Difference (%)',xlabel='Scalar (ell,m)')
    err.set_xticks(x,[f"({r['ell']},{r['m']})" for r in rows]);err.grid(axis='y',alpha=.2)
    fig.savefig(stem.with_suffix('.png'),dpi=170);plt.close(fig)
    print(json.dumps({k:report[k] for k in ('finite_infinity_total','finite_horizon_total','sum_figure6_ell_le6')}))


if __name__=='__main__':main()
