"""Compare completed baseline/dense matching and selected actual sources."""
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    root=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    prefix='metric_tensor_matching_q24_smooth_m1_r20_a0.877153027595'
    paths=[root/(prefix+suffix+'.json') for suffix in ('','_denseangular')]
    paths.append(root/'dense_source_mg2_selected_nodes.json')
    old,new,source=[json.loads(p.read_text()) for p in paths]
    for field in ('jet_order','smooth_polar_tests','r0','a','m','quadrature','epsilon'):
        if old[field]!=new[field]:raise ValueError(f'Matching parameters differ: {field}')
    if old['cases'][-1]['ellmax']!=18 or new['cases'][-1]['ellmax']!=18:
        raise ValueError('Both matching calculations must reach ellmax=18')
    if not source['all_requested_nodes_completed']:raise ValueError('Incomplete selected-source run')
    a,b=old['cases'][-1]['limits'][-1],new['cases'][-1]['limits'][-1]
    checks={key:dict(baseline=a[key],dense=b[key],reduction_factor=a[key]/b[key])
            for key in ('maximum_value_jump','maximum_derivative_error')}
    fig,axes=plt.subplots(1,2,figsize=(10,4.6),layout='constrained')
    for data,label,color in ((old,'Baseline angular solver','.45'),(new,'Dense angular solver','#1764ab')):
        ell=[case['ellmax'] for case in data['cases']]
        for ax,key in zip(axes,checks):
            ax.semilogy(ell,[case['limits'][-1][key] for case in data['cases']],
                        'o-',color=color,label=label,ms=3,lw=1)
    for ax,title in zip(axes,('Metric value jump','Derivative jump error')):
        ax.set(xlabel=r'Metric $\ell_{\max}$',ylabel='Maximum absolute component',title=title,xticks=range(2,19,2))
        ax.grid(alpha=.2)
    axes[0].legend(fontsize=8)
    fig.suptitle('Kerr point-source matching: r_p=20M, m_g=1\nResiduals remain; this is not a flux-error bound')
    stem=root/'dense_angular_validation'
    fig.savefig(stem.with_suffix('.png'),dpi=160);plt.close(fig)
    result=dict(status='completed_diagnostics_not_flux_convergence',
        inputs_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        matching_at_L18=checks,selected_source_comparisons=source['completed'],
        limitations=['Dense backend does not eliminate the high-L matching error floor',
                     'Selected source nodes do not determine a corrected integrated flux',
                     'Default production angular backend remains unchanged'])
    stem.with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(checks,indent=2),flush=True)


if __name__=='__main__':main()
