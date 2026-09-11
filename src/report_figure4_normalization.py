"""Cross-check Fig.4 prefactors against Fig.2 and direct GW modes."""
import hashlib
import json
from pathlib import Path
from report_gravitational_flux import mode_flux


def main():
    folder=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    names=['paper_figure4_selected_radii.json','paper_figure2_selected_radii.json',
           'gravitational_flux_r20_L12.json',
           'flux_coverage_L18_nt18_i6_h5_f5_m0go4000_m2go4000.json']
    f4,f2,gw,scalar=[json.loads((folder/n).read_text()) for n in names]
    alpha=.3;q=1e-6;eta=.1;a=gw['parameters']['a']
    f2rows={r['rp']:r for r in f2['rows']};rows=[]
    factors=dict(prose_eta=q*q*eta*alpha**6,caption_epsilon=q*q*.1*alpha**3,
                 extra_alpha6_hypothesis=q*q*eta*alpha**12)
    for row in f4['rows']:
        rp=row['rp'];direct=mode_flux(a,rp,2,2);scales={}
        for boundary in ('infinity','horizon_magnitude'):
            if rp in f2rows:
                source=f2rows[rp]['flux'][boundary]['plotted_flux']
                measured=row['scalar'][boundary]['value']
                scales[boundary]=dict(figure4_over_figure2=measured/source,
                    measured_over_predictions={k:measured/(source*v) for k,v in factors.items()})
        rows.append(dict(rp=rp,scalar_scalings=scales,
            inset_point_visible={b:row['ratio'][b]['inside_axes'] for b in ('infinity','horizon_magnitude')},
            single_positive_22_gw_per_q2=direct,
            figure4_gw_over_single_22={b:row['gravitational'][b]['value']/(q*q*abs(direct['infinity' if b=='infinity' else 'horizon']))
                for b in ('infinity','horizon_magnitude')},
            inset_over_main_flux_ratio={b:row['ratio'][b]['value']/(row['scalar'][b]['value']/row['gravitational'][b]['value'])
                for b in ('infinity','horizon_magnitude')}))
    if any(scalar[b]['finite_resolution_total'] is None for b in ('infinity','horizon')):
        raise ValueError('Require complete finite scalar totals')
    numerical={b:dict(scalar_physical=q*q*eta*scalar[b]['finite_resolution_total'],
        gravitational_physical=q*q*gw['totals'][b],
        scalar_over_gravitational=eta*scalar[b]['finite_resolution_total']/gw['totals'][b])
        for b in ('infinity','horizon')}
    result=dict(status='figure4_cross_figure_normalization_audit_not_reproduction',
        input_sha256={n:hashlib.sha256((folder/n).read_bytes()).hexdigest() for n in names},
        parameters=dict(alpha=alpha,q=q,eta=eta,a=a),predetermined_factors=factors,
        reference_rows=rows,r20_numerical=numerical,
        limitations=['Figure values are digitized and have no rigorous error bound',
            'Agreement with a single GW mode is evidence about the plotted curve, not proof of author code',
            'Extra alpha^6 is a tested algebraic hypothesis, not a correction applied to any solver',
            'Scalar source and metric convergence remain unresolved; this is not a physical flux prediction'])
    path=folder/'figure4_normalization_audit.json';path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
