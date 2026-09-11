"""Reproject the completed finite even spheroidal field onto spherical modes.

Includes the scalar (0,0) field: although excluded from wake plots, its
spheroidal angular function has small nonzero higher spherical components.
"""
import hashlib
import json
from pathlib import Path
import numpy as np
from environment_source import angular_mode


def main():
    folder=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    coverage_path=folder/'flux_coverage_L18_nt18_i12_h12_f12_m0go4000_m2go4000.json'
    coverage=json.loads(coverage_path.read_text());inputs=[];common=None
    required={(ell,m) for ell in range(0,13,2) for m in range(-ell,ell+1,2)}
    records={(row['ell'],row['m']):row for row in coverage['horizon']['modes']}
    for ell,m in sorted(required):
        row=records[(ell,m)]
        if 'flux' not in row:raise ValueError('Complete even input shells including ell=0 are required')
        path=folder/row['file'];data=json.loads(path.read_text());p=data['parameters'];g=p['metric']
        if (data['status']!='truncated_single_mode_not_converged'
            or data.get('flux_validity')=='historical_unreliable_boundary_result'):
            raise ValueError('Incomplete or invalid source')
        if (p['scalar_ell'],p['scalar_m'])!=(ell,m):raise ValueError('Inventory mode mismatch')
        physics=(p['alpha'],g['a'],g['orbital_radius'],p['cloud_mass'])
        if common is not None and physics!=common:raise ValueError('Incompatible input physics')
        common=physics
        point=[r for r in data['radial_response_at_panel_boundaries'] if r['r']==g['orbital_radius']]
        if len(point)!=1:raise ValueError('Missing exact particle radial value')
        inputs.append(dict(ell=ell,m=m,omega=p['omega'],radial=point[0]['field'],
                           file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    alpha,a,r0,mass=common
    if mass!=1:raise ValueError('Unit cloud mass required for this diagnostic')
    shells_by_order={};reconstruction={}
    for order in (80,120):
        x,w=np.polynomial.legendre.leggauss(order);theta=np.arccos(x)
        spherical={};point_values={};original=0j
        shells={L:0j for L in range(0,33,2)}
        for row in inputs:
            ell,m,omega=row['ell'],row['m'],row['omega']
            source=angular_mode(theta,ell,m,a*a*(omega*omega-alpha*alpha))[0]
            np.testing.assert_allclose(2*np.pi*np.dot(w,source*source),1.,rtol=1e-12)
            radial=complex(*row['radial'])
            original+=radial*angular_mode(np.pi/2,ell,m,a*a*(omega*omega-alpha*alpha))[0]
            for L in range(abs(m),33,2):
                key=(L,m)
                if key not in spherical:
                    size=max(20,L-abs(m)+1)
                    spherical[key]=angular_mode(theta,L,m,0.,size=size)[0]
                    point_values[key]=angular_mode(np.pi/2,L,m,0.,size=size)[0]
                b=2*np.pi*np.dot(w,spherical[key]*source)
                shells[L]+=radial*b*point_values[key]
        reconstructed=sum(shells.values())
        np.testing.assert_allclose(reconstructed,original,rtol=1e-10,atol=1e-14)
        reconstruction[order]=dict(original=[original.real,original.imag],
            projected=[reconstructed.real,reconstructed.imag],relative_difference=abs(reconstructed-original)/abs(original))
        shells_by_order[order]=shells
    rows=[]
    for L in range(0,33,2):
        value=shells_by_order[120][L]
        rows.append(dict(ell=L,spherical_field=[value.real,value.imag],
            absolute_quadrature_difference=abs(value-shells_by_order[80][L]),
            scaled_local_coefficient=[-(L+.5)**2*value.real/alpha**3,-(L+.5)**2*value.imag/alpha**3]))
    result=dict(status='finite_spheroidal_field_reprojection_not_infinite_mode_convergence',
        coverage_sha256=hashlib.sha256(coverage_path.read_bytes()).hexdigest(),
        parameters=dict(alpha=alpha,a=a,rp=r0,cloud_mass=mass,input_spheroidal_ellmax=12,
            output_spherical_ellmax=32,quadrature_orders=[80,120],theta=np.pi/2,phi=0,time=0),
        inputs=inputs,rows=rows,reconstruction=reconstruction,
        limitations=['All even spheroidal inputs 0..12 included, no static odd-m field needed',
            'Spherical ell>12 output is angular mixing of retained inputs, not new physical solved modes',
            'Missing spheroidal ell>12 can contribute to every output shell',
            'No infinite-mode or source-grid convergence claim'])
    (folder/'particle_spherical_projection.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(reconstruction=reconstruction,high_modes=[r for r in rows if r['ell'] in (6,8,10,12)])))


if __name__=='__main__':main()
