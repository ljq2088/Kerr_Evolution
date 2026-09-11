"""Compare point-source matching diagnostics without inferring flux errors."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    root=Path(__file__).resolve().parents[1]
    directory=root/'docs/environment_reproduction'
    prefix='metric_tensor_matching_q{}_smooth_m1_r41.6_a0.877153027595{}.json'
    runs=[]
    fig,axes=plt.subplots(1,2,figsize=(10,4.2),constrained_layout=True)
    for n,suffix,label in ((12,'','12 nodes, double'),(16,'','16 nodes, double'),
                           (12,'_extended','12 nodes, extended jets')):
        path=directory/prefix.format(n,suffix)
        data=json.loads(path.read_text())
        scale=float(np.max(abs(np.asarray(data['expected_derivative_jump']))))
        cases=[c for c in data['cases'] if c['ellmax']<=8]
        ell=[c['ellmax'] for c in cases]
        value=[c['limits'][-1]['maximum_value_jump'] for c in cases]
        derivative=[c['limits'][-1]['maximum_derivative_error'] for c in cases]
        runs.append(dict(label=label,source=path.name,ellmax=ell,
            maximum_value_jump=value,maximum_derivative_error=derivative,
            derivative_target_maximum=scale,
            final_derivative_error_over_maximum_target=derivative[-1]/scale))
        for ax,values in zip(axes,(value,derivative)):
            ax.semilogy(ell,values,'o-',label=label,markersize=4)
    for ax,title in zip(axes,('Metric value jump','Derivative jump error')):
        ax.set(xlabel=r'$\ell_{g,\max}$',ylabel='Maximum absolute component',title=title)
        ax.grid(alpha=.25)
    axes[0].legend(fontsize=8)
    fig.suptitle(r'Point-source diagnostic: $r_p=41.6M$, $m_g=1$')
    output=root/'outputs/environment_validation'
    output.mkdir(parents=True,exist_ok=True)
    fig.savefig(output/'matching_precision.png',dpi=150)
    plt.close(fig)
    (directory/'matching_precision_comparison.json').write_text(json.dumps(dict(
        status='diagnostic_not_a_flux_error_bound',runs=runs,
        conclusion='Increasing quadrature and jet precision alone does not remove the observed error floor; radial and harmonic inputs remain double precision.'),indent=2)+'\n')


if __name__=='__main__':
    main()
