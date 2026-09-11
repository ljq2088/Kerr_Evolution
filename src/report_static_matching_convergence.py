"""Summarize static junction evidence without claiming a physical-boundary solution."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    root=Path(__file__).resolve().parents[1]
    folder=root/'docs/environment_reproduction'
    names=['static_tetrad_a0.6_r6_L8_q20_j8_free8_paper.json',
           'static_tetrad_a0.6_r6_L12_q20_j12_free12_paper.json',
           'static_tetrad_a0.6_r6_L12_q20_j12_free12_paper_eps5e-06.json']
    cases=[json.loads((folder/name).read_text()) for name in names]
    for case in cases:
        assert case['sample_metadata']['a']==.6 and case['sample_metadata']['r0']==6.
        assert not case['derivative_conditions_used_in_fit']
        assert case['charge_conditions']=='exact_elimination'
    coefficients=dict(zip(cases[-1]['basis_labels'],cases[-1]['coefficients']))
    paper=dict(B=-.17354798,D=.19491814,F=.025072823)
    even=(2,4,6,8)
    ratio={ell:cases[1]['scaled_derivative_residual_by_degree'][ell]/cases[2]['scaled_derivative_residual_by_degree'][ell] for ell in even}
    result=dict(status='local_low_multipole_matching_verified_boundaries_and_full_sum_incomplete',
        inputs=names,reference='https://arxiv.org/abs/2306.16459v3',
        reference_table='Jumps in radial functions and completion mode coefficients for a=0.6M, r0=6M',
        independent_completion_comparison={k:dict(computed=coefficients[k],paper_rounded=v,absolute_difference=abs(coefficients[k]-v)) for k,v in paper.items()},
        not_compared='C and free scalar jumps depend on the particular/homogeneous gauge convention; E/G are imposed exactly, not independent validation',
        derivative_error_ratio_epsilon_over_tenth=ratio,
        verified_test_degrees=list(range(9)),
        low_degree_maximum_value=float(max(cases[2]['scaled_value_residual_by_degree'][:9])),
        low_degree_maximum_derivative=float(max(cases[2]['scaled_derivative_residual_by_degree'][:9])))
    (folder/'static_matching_convergence.json').write_text(json.dumps(result,indent=2)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    labels=[r'$L=8,\ \epsilon=5\times10^{-5}$',r'$L=12,\ \epsilon=5\times10^{-5}$',r'$L=12,\ \epsilon=5\times10^{-6}$']
    for case,label in zip(cases,labels):
        for ax,key in zip(axes,('scaled_value_residual_by_degree','scaled_derivative_residual_by_degree')):
            values=np.asarray(case[key])
            ax.semilogy(np.arange(len(values)),np.maximum(values,1e-17),'.-',label=label)
    for ax,title in zip(axes,('All five tetrad components: continuity','Independent radial derivative jumps')):
        ax.set(xlabel='Angular test degree j',ylabel='Maximum scaled residual',title=title)
        ax.grid(alpha=.25);ax.set_xticks(range(0,15,2));ax.axvspan(9.5,14,alpha=.08,color='red')
    axes[0].legend(fontsize=8)
    fig.suptitle('Static Kerr local matching: a/M=0.6, r0/M=6\nShaded high-degree region is not converged; physical boundaries are not imposed',fontsize=11)
    output=root/'outputs/environment_validation'
    output.mkdir(exist_ok=True,parents=True)
    fig.savefig(output/'static_matching_convergence.png',dpi=170)
    fig.savefig(output/'static_matching_convergence.pdf')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
