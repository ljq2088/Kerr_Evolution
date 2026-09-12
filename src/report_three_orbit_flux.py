"""Compare complete finite alpha=.3 Kerr sums with digitized Figure 2 values."""
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    folder=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    reference_path=folder/'paper_figure2_selected_radii.json'
    reference=json.loads(reference_path.read_text())
    refs={row['rp']:row for row in reference['rows']}
    names={10:'flux_coverage_L18_nt18_i6_h5_f5_rp10.json',
           20:'flux_coverage_L18_nt18_i6_h5_f5_m0go4000_m2go4000.json',
           30:'flux_coverage_L18_nt18_i6_h5_f5_rp30.json'}
    rows=[];inputs={reference_path.name:hashlib.sha256(reference_path.read_bytes()).hexdigest()}
    for radius,name in names.items():
        path=folder/name;coverage=json.loads(path.read_text());p=coverage['parameters']
        if p['alpha']!=.3 or p['rp']!=radius or p['metric_ellmax']!=18:
            raise ValueError('Wrong coverage physics or metric truncation')
        result=dict(rp=radius,fluxes={})
        inputs[name]=hashlib.sha256(path.read_bytes()).hexdigest()
        for boundary,count in (('infinity',9),('horizon',18)):
            data=coverage[boundary]
            if data['required_count']!=count or data['computed_count']!=count:
                raise ValueError('Finite range incomplete')
            values=[]
            for channel in data['modes']:
                source=folder/channel['file'];actual=json.loads(source.read_text())
                if actual['status']!='truncated_single_mode_not_converged':raise ValueError('Incomplete source')
                if actual['parameters'].get('background') is not None:raise ValueError('Wrong background')
                value=actual['flux'][boundary]['orbital_energy']
                if value!=channel['flux'][boundary]['orbital_energy']:raise ValueError('Stale coverage flux')
                values.append(value);inputs[source.name]=hashlib.sha256(source.read_bytes()).hexdigest()
            total=sum(values)
            if total!=data['finite_resolution_total']:raise ValueError('Coverage sum mismatch')
            key='infinity' if boundary=='infinity' else 'horizon_magnitude'
            target=refs[radius]['flux'][key]['per_q2_cloud_mass']
            result['fluxes'][boundary]=dict(computed_signed=total,paper_magnitude=target,
                relative_magnitude_difference=abs(total)/target-1,channels=count)
        rows.append(result)
    report=dict(status='finite_three_orbit_comparison_not_paper_reproduction',alpha=.3,
        normalization='per q^2 (Mc/M)',rows=rows,inputs=inputs,amplitude_fitted=False,
        limitations='Only three sampled orbits. Paper values are digitized without a rigorous readout error bound. Horizon reference gives magnitude only. Source, metric and all-mode convergence remain unresolved.')
    stem=folder/'figure2_three_orbit_comparison'
    Path(str(stem)+'.json').write_text(json.dumps(report,indent=2)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(10,4.5),layout='constrained')
    x=[row['rp'] for row in rows]
    for ax,boundary in zip(axes,('infinity','horizon')):
        computed=[abs(row['fluxes'][boundary]['computed_signed']) for row in rows]
        paper=[row['fluxes'][boundary]['paper_magnitude'] for row in rows]
        ax.semilogy(x,computed,'o',label='Computed finite mode sum')
        ax.semilogy(x,paper,'s',fillstyle='none',label='Paper Figure 2 readout')
        for radius,y,row in zip(x,computed,rows):
            delta=row['fluxes'][boundary]['relative_magnitude_difference']*100
            ax.annotate(f'{delta:+.2f}%',(radius,y),xytext=(0,8),textcoords='offset points',ha='center')
        ax.set(xlabel=r'$r_p/M$',ylabel=r'$|F|/[q^2(M_c/M)]$',
               title='Infinity' if boundary=='infinity' else 'Horizon magnitude',xticks=x,xlim=(7,33))
        ax.grid(alpha=.25);ax.margins(y=.2)
    axes[0].legend(fontsize=8)
    fig.suptitle('Kerr, alpha=0.3: three sampled orbits (not a converged radial scan)')
    fig.savefig(Path(str(stem)+'.png'),dpi=160);plt.close(fig)
    print(json.dumps(rows),flush=True)


if __name__=='__main__':main()
