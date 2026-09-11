"""Show calculated and missing Figure 6 channels without a partial total."""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('coverage',type=Path)
    args=parser.parse_args()
    root=args.coverage.parent
    coverage=json.loads(args.coverage.read_text())
    reference_path=root/'paper_v1_flux_markers.json'
    reference=json.loads(reference_path.read_text())
    targets={(r['ell'],r['m']):r['plotted_flux'] for r in reference['markers']}
    section=coverage['infinity'];computed=[];missing=[]
    for row in section['modes']:
        if 'flux' not in row:
            missing.append([row['ell'],row['m']]);continue
        value=row['flux']['infinity']['orbital_energy']
        if value<=0:raise ValueError('Expected a positive propagating orbital flux')
        target=targets.get((row['ell'],row['m']))
        computed.append(dict(ell=row['ell'],m=row['m'],file=row['file'],computed=value,
                             digitized=target,ratio=value/target if target is not None else None))
    count=len(computed);required=section['required_count']
    fig,(ax,ratio)=plt.subplots(2,1,figsize=(9,7),sharex=True,
                              gridspec_kw={'height_ratios':[2,1]},layout='constrained')
    ax.semilogy([r['ell'] for r in reference['markers']],
                [r['plotted_flux'] for r in reference['markers']],
                'o',ms=5,mfc='none',mec='.5',label='Paper v1: digitized markers')
    colors=plt.get_cmap('tab20')
    for index,m in enumerate(sorted({r['m'] for r in computed})):
        rows=sorted((r for r in computed if r['m']==m),key=lambda r:r['ell'])
        ax.semilogy([r['ell'] for r in rows],[r['computed'] for r in rows],
                    'x-',lw=.8,ms=6,color=colors(index%20),label=f'Computed m={m}')
        matched=[r for r in rows if r['ratio'] is not None]
        ratio.plot([r['ell'] for r in matched],[100*(r['ratio']-1) for r in matched],
                   'o-',ms=4,lw=.8,color=colors(index%20))
    ax.set(ylim=(1e-17,1e-4),ylabel=r'Infinity orbital flux per $q^2(M_c/M)$')
    ax.legend(ncol=2,fontsize=8);ax.grid(alpha=.2)
    ratio.axhline(0,color='.4',lw=.8)
    ratio.set(xlabel=r'Scalar $\ell$',ylabel='Difference (%)',xticks=range(2,13))
    ratio.grid(alpha=.2)
    fig.suptitle(f'Figure 6 comparison: {count}/{required} channels computed\n'
                 'Kerr, alpha=0.3, r_p=20M; finite resolution, not converged')
    stem=root/f'figure6_progress_L{coverage["parameters"]["metric_ellmax"]}'
    fig.savefig(stem.with_suffix('.png'),dpi=170);plt.close(fig)
    result=dict(status='partial_figure6_comparison' if missing else 'finite_figure6_comparison_not_converged',
                parameters=coverage['parameters'],computed_count=count,required_count=required,
                modes=computed,missing_modes=missing,finite_total=section['finite_resolution_total'],
                inputs_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (args.coverage,reference_path)},
                limitations=['No missing mode is assigned zero; no partial sum is called a total',
                             'Values outside displayed flux range remain in this JSON',
                             'Paper markers are digitized, not author numerical data',
                             'Normalization inferred from Figures 2 and 6 closure',
                             'Angular repeatability and source-grid errors remain under investigation'])
    stem.with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(stem.with_suffix('.png'),count,required,flush=True)


if __name__=='__main__':main()
