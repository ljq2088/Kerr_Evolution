"""Compare complete r20 finite totals against original Fig.2 vector curves."""
import hashlib
import json
from pathlib import Path
import numpy as np


def main():
    folder=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    rows=[];inputs=[];ratios=[]
    for background,suffix in [('schwarzschild','_schwarzschild_frozen'),('kerr','')]:
        cp=folder/f'flux_coverage_L18_nt18_i6_h5_f5{suffix}.json'
        rp=folder/f'figure2_{background}_r20_reference.json'
        data=json.loads(cp.read_text());ref=json.loads(rp.read_text());p=data['parameters']
        if p['alpha']!=.3 or p['rp']!=20 or p['metric_ellmax']!=18 or p['nt']!=18:
            raise ValueError('Wrong comparison parameters')
        if ref['background']!=background or ref['rows'][0]['rp']!=20:
            raise ValueError('Wrong reference curve')
        checked={}
        for boundary,key,count in [('infinity','infinity',9),('horizon','horizon_magnitude',18)]:
            block=data[boundary]
            if block['required_count']!=count or block['computed_count']!=count or block['finite_resolution_total'] is None:
                raise ValueError('Incomplete finite flux inventory')
            modes=set();total=0.
            for row in block['modes']:
                mode=(row['ell'],row['m'])
                if mode in modes:raise ValueError('Duplicate mode')
                modes.add(mode)
                path=cp.parent/row['file'];response=json.loads(path.read_text())
                q=response['parameters']
                if response['status']!='truncated_single_mode_not_converged' or response['flux']!=row['flux']:
                    raise ValueError('Response differs from inventory')
                if (q['scalar_ell'],q['scalar_m'])!=mode:
                    raise ValueError('Response mode differs')
                if not np.isfinite([x for s in response['samples'] for x in [s['r'],s['weight'],*s['source']]]).all():
                    raise ValueError('Nonfinite source')
                value=response['flux'][boundary]['orbital_energy']
                if not np.isfinite(value):raise ValueError('Nonfinite flux')
                total+=value;checked[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
            np.testing.assert_allclose(total,block['finite_resolution_total'],rtol=1e-14,atol=0)
            reference=ref['rows'][0]['flux'][key];paper=reference['per_q2_cloud_mass']
            x=ref['rows'][0]['pdf_x'];(x0,y0),(x1,y1)=reference['bracketing_pdf_points']
            dy=(ref['ygrid'][-1]-ref['ygrid'][0])/3
            f0,f1=[10**(-1-(y-ref['ygrid'][0])/dy) for y in (y0,y1)]
            linear_flux=f0+(f1-f0)*(x-x0)/(x1-x0)
            rows.append(dict(background=background,boundary=boundary,
                signed_computed_per_q2_eta=total,paper_magnitude_per_q2_eta=paper,
                computed_magnitude_per_q2_epsilon2=abs(total)/.3**6,
                paper_magnitude_per_q2_epsilon2=reference['plotted_flux'],
                relative_magnitude_difference_percent=100*(abs(total)/paper-1),
                interpolation_method_sensitivity_percent=100*(linear_flux/reference['plotted_flux']-1)))
        selected=rows[-2:]
        actual_ratio=abs(selected[1]['signed_computed_per_q2_eta']/selected[0]['signed_computed_per_q2_eta'])
        paper_ratio=selected[1]['paper_magnitude_per_q2_eta']/selected[0]['paper_magnitude_per_q2_eta']
        ratios.append(dict(background=background,computed_horizon_over_infinity=actual_ratio,
            paper_horizon_over_infinity=paper_ratio,relative_difference_percent=100*(actual_ratio/paper_ratio-1)))
        inputs.append(dict(background=background,coverage=cp.name,reference=rp.name,
            coverage_sha256=hashlib.sha256(cp.read_bytes()).hexdigest(),reference_sha256=hashlib.sha256(rp.read_bytes()).hexdigest(),
            response_sha256=checked))
    result=dict(status='finite_totals_vs_digitized_paper_not_convergence',alpha=.3,rp=20,
        rows=rows,normalization_invariant_ratios=ratios,inputs=inputs,
        limitations=['Reference values come from original PDF vector interpolation, not author tables.',
            'Interpolation method sensitivity is not a bound on digitization or paper discretization error.',
            'Horizon comparison uses magnitudes; signed computed horizon fluxes remain negative.',
            'Single radius, finite modal truncation, and frozen Schwarzschild decay do not establish full convergence.'])
    output=folder/'background_r20_paper_comparison.json';output.write_text(json.dumps(result,indent=2)+'\n')
    for row in rows:print(row)
    print(ratios)


if __name__=='__main__':main()
