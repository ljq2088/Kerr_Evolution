"""Reintegrate saved actual-source data and audit a horizon-tail correction."""
import argparse
import json
from pathlib import Path
import numpy as np
from environment_radial import RadialGreen
from environment_cloud import mode_flux
from environment_horizon_tail import horizon_tail


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('inputs',nargs='+',type=Path)
    args=parser.parse_args()
    cases=[]
    for path in args.inputs:
        data=json.loads(path.read_text());p=data['parameters'];m=p['metric']
        if 'flux' not in data:raise ValueError('Source integration must finish first')
        if p.get('background') is not None or p['scalar_m']-m['m_g']!=1:
            raise ValueError('This phase prescription requires the stationary threshold |211> cloud')
        if data.get('flux_validity')=='historical_unreliable_boundary_result':
            raise ValueError('Recompute the radial boundary before using this response')
        a,mu,omega=m['a'],p['alpha'],p['omega']
        green=RadialGreen(a,mu,omega,p['scalar_ell'],p['scalar_m'],rmax=p['green_outer_radius'],
                          offset=p['green_horizon_offset'],rtol=1e-11,infinity_method=p.get('infinity_method','series'))
        r=np.array([s['r'] for s in data['samples']]);J=np.array([complex(*s['source']) for s in data['samples']])
        weights=np.array([s['weight'] for s in data['samples']])
        integrand=green.upsol.sol(r)[0]*J/green.w0
        base=np.dot(weights,integrand)
        np.testing.assert_allclose(base,complex(*data['z_h']),rtol=3e-10,atol=2e-13)
        cutoff=p['source_panels'][0]-green.rp
        gamma=(2*green.rp*omega-a*p['scalar_m'])/(2*np.sqrt(1-a*a))
        omega_c=omega-m['m_g']/(m['orbital_radius']**1.5+a)
        variations=[]
        for order in (0,1,2):
            for window in (5.,10.,20.):
                if cutoff*window>=min(2*np.sqrt(1-a*a),m['orbital_radius']-green.rp):
                    variations.append(dict(order=order,window=window,
                        unavailable='Window reaches the inner-horizon distance or the particle-source radius'))
                    continue
                try:correction,report=horizon_tail(r,integrand,green.rp,gamma,cutoff,order=order,window=window)
                except ValueError as error:
                    variations.append(dict(order=order,window=window,unavailable=str(error)));continue
                corrected=base+correction
                flux=mode_flux(omega,p['scalar_m'],omega_c,1,mu,a,0j,corrected)['horizon']['orbital_energy']
                variations.append(dict(**report,correction=[correction.real,correction.imag],
                    corrected_z_h=[corrected.real,corrected.imag],corrected_horizon_flux=flux))
        cases.append(dict(input=path.name,cutoff=cutoff,horizon_order=p.get('horizon_quadrature_order'),
            base_horizon_flux=data['flux']['horizon']['orbital_energy'],gamma=gamma,variations=variations))
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/horizon_tail_audit.json'
    out.write_text(json.dumps(dict(status='horizon_tail_extrapolation_diagnostic_not_converged',cases=cases),indent=2)+'\n')
    for case in cases:
        print(case['input'],flush=True)
        for v in case['variations']:
            if 'corrected_horizon_flux' in v:
                print(v['order'],v['window'],v['corrected_horizon_flux'],v['held_out_integrand_residual'],flush=True)


if __name__=='__main__':main()
