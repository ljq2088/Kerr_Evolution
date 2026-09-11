"""Coherent field sums at the particle; incomplete ell shells are never summed."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from environment_source import angular_mode


def encode(z):return [float(z.real),float(z.imag)]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('coverage',type=Path)
    args=parser.parse_args()
    folder=args.coverage.parent
    coverage=json.loads(args.coverage.read_text())
    rows=coverage['field']['modes'];multipoles=[];inputs={};common=None
    for ell in sorted({r['ell'] for r in rows}):
        required=list(range(-ell,ell+1,2));shell={r['m']:r for r in rows if r['ell']==ell}
        missing=[m for m in required if m not in shell or 'flux' not in shell[m]]
        if missing:
            multipoles.append(dict(ell=ell,missing_m=missing,complete=False));continue
        components=[]
        for m in required:
            path=folder/shell[m]['file'];data=json.loads(path.read_text());p=data['parameters'];g=p['metric']
            if data['status']!='truncated_single_mode_not_converged' or p['scalar_ell']!=ell or p['scalar_m']!=m:
                raise ValueError('Inventory/source mismatch')
            if data.get('flux_validity')=='historical_unreliable_boundary_result':raise ValueError('Invalid boundary')
            if p['cloud_mass']!=1:raise ValueError('This diagnostic requires unit cloud mass')
            physical=(p['alpha'],g['a'],g['orbital_radius'],p['cloud_mass'],data['flux_scaling'])
            if common is not None and physical!=common:raise ValueError('Incompatible physical normalization')
            common=physical
            points=[row for row in data['radial_response_at_panel_boundaries'] if row['r']==g['orbital_radius']]
            if len(points)!=1:raise ValueError('Need one exact orbit panel boundary')
            radial=complex(*points[0]['field'])
            angular=angular_mode(np.pi/2,ell,m,g['a']**2*(p['omega']**2-p['alpha']**2))[0]
            components.append(dict(m=m,radial=encode(radial),equatorial_angular_value=float(angular),
                particle_contribution=encode(radial*angular),file=path.name))
            inputs[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        radial_sum=sum(complex(*row['radial']) for row in components)
        particle_sum=sum(complex(*row['particle_contribution']) for row in components)
        multipoles.append(dict(ell=ell,complete=True,missing_m=[],components=components,
            radial_coefficient_sum=encode(radial_sum),particle_field_sum=encode(particle_sum),
            abs_radial_coefficient_sum=abs(radial_sum),abs_particle_field_sum=abs(particle_sum),
            abs_radial_sum_per_epsilon_q=abs(radial_sum)/common[0]**3,
            abs_particle_sum_per_epsilon_q=abs(particle_sum)/common[0]**3))
    result=dict(status='finite_particle_field_not_figure7_reproduction',
        evaluation=dict(r=coverage['parameters']['rp'],theta=np.pi/2,phi=0,time=0),
        coverage_sha256=hashlib.sha256(args.coverage.read_bytes()).hexdigest(),inputs_sha256=inputs,
        multipoles=multipoles,complete_shells=[r['ell'] for r in multipoles if r['complete']],
        limitations=['Original finite Gauss source integrals at exact panel boundary; no interpolation to particle',
            'Complex amplitudes summed coherently before absolute value',
            'Radial-only sum and angular-evaluated physical field are distinct quantities',
            'alpha^-3 conversion uses epsilon=alpha^3 sqrt(Mc/M); not a fitted scale',
            'No missing m channel is set to zero; no asymptotic slope claim from incomplete sequence',
            'Metric truncation basis and source convergence remain unverified'])
    path=folder/f"particle_field_L{coverage['parameters']['metric_ellmax']}.json"
    path.write_text(json.dumps(result,indent=2)+'\n')
    markers_path=folder/'paper_figure7_markers.json'
    if markers_path.exists():
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        markers=json.loads(markers_path.read_text())['markers']
        complete=[row for row in multipoles if row['complete']]
        fig,ax=plt.subplots(figsize=(8,5),layout='constrained')
        ax.loglog([row['ell'] for row in markers],[row['plotted_value'] for row in markers],
                  'D',color='.45',label='Paper Fig.7 vector markers')
        ax.loglog([row['ell'] for row in complete],[row['abs_particle_sum_per_epsilon_q'] for row in complete],
                  'o',color='#1667b1',label='Angular-evaluated field / (epsilon q)')
        ax.loglog([row['ell'] for row in complete],[row['abs_radial_sum_per_epsilon_q'] for row in complete],
                  'x',color='#bc6537',label='Radial-coefficient sum / (epsilon q)')
        ax.set(xlabel='Scalar ell',ylabel='Absolute coherent amplitude',
               title='Particle at r=20M: definition audit, not a Fig.7 reproduction')
        ax.set_xticks(range(2,13),labels=[str(i) for i in range(2,13)])
        ax.legend();ax.grid(alpha=.2)
        fig.savefig(path.with_suffix('.png'),dpi=160);plt.close(fig)
    print(json.dumps({'complete_shells':result['complete_shells']}))


if __name__=='__main__':main()
