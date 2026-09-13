"""Isolated coherent scalar00 flux sensitivity to dipole kappa boundary data.

Other metric sectors, source quadrature and scalar Green function stay fixed.
No metric cache entries are modified by this diagnostic.
"""
import hashlib
import json
from pathlib import Path
import numpy as np
from environment_source import ThresholdCloud, angular_mode, connection
from environment_radial import RadialGreen
from environment_trace_variation import TraceMassVariation
from lorenz_kappa import kappa_jet
from lorenz_ghp import KerrGHP


def scalar_hessian(k, omega, m, gamma):
    f=k.value
    grad=np.array([-1j*omega*f,k.derivative(0).value,k.derivative(1).value,1j*m*f])
    partial=np.empty((4,4),complex)
    partial[0]=-1j*omega*grad;partial[:,0]=partial[0]
    partial[3]=1j*m*grad;partial[:,3]=partial[3]
    partial[1,1]=k.derivative(0).derivative(0).value
    partial[2,2]=k.derivative(1).derivative(1).value
    partial[1,2]=partial[2,1]=k.derivative(0).derivative(1).value
    return partial-np.einsum('cab,c->ab',gamma,grad)


def main():
    folder=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    path=folder/'forced_mode_nr8_nt18_L1_mg-1_sl0_inner0.0005_outer320_log_h32.json'
    raw=path.read_bytes();data=json.loads(raw);p=data['parameters']
    assert data['status']=='truncated_single_mode_not_converged'
    assert p['scalar_ell']==p['scalar_m']==0 and p['metric']['ellmax']==1
    assert p['metric']['m_g']==-1
    cloud=ThresholdCloud(alpha=p['alpha']);r0=p['metric']['orbital_radius']
    assert abs(cloud.a-p['metric']['a'])<1e-13
    omega_g=1/(r0**1.5+cloud.a)
    green=RadialGreen(cloud.a,p['alpha'],p['omega'],0,0,
        rmax=p['green_outer_radius'],offset=p['green_horizon_offset'],rtol=1e-11,
        infinity_method=p.get('infinity_method','series'))
    radii=np.array([s['r'] for s in data['samples']]);weights=np.array([s['weight'] for s in data['samples']])
    source=np.array([complex(*s['source']) for s in data['samples']])
    kernel=weights*green.upsol.sol(radii)[0]/green.w0
    saved=complex(*data['z_h']);recomputed=kernel@source
    if abs(recomputed/saved-1)>1e-9:raise ValueError('Saved amplitude not reproduced')
    x,w=np.polynomial.legendre.leggauss(p['angular_order']);theta=np.arccos(x)
    angular=angular_mode(theta,0,0,cloud.a**2*(p['omega']**2-cloud.mu**2))[0]
    configurations=[2000.,4000.,8000.]
    variations=[TraceMassVariation(r0,cloud.a,1,1,rmax=outer) for outer in configurations]
    corrections=np.zeros((len(configurations),len(radii)),complex)
    out=folder/'dipole_kappa_flux_boundary_audit.json'
    for i,r in enumerate(radii):
        for t,xx,ww,ss in zip(theta,x,w,angular):
            g=KerrGHP(float(r),float(t),cloud.a,omega=omega_g,m=1,order=2)
            baseline=kappa_jet(g,r0,1)
            inv,gamma=connection(r,t,cloud.a)
            cloud_hessian=cloud.hessian(r,t)[1]
            raised=inv@cloud_hessian@inv
            factor=2*np.pi*ww*ss*(r*r+cloud.a**2*xx*xx)
            for j,variation in enumerate(variations):
                delta=variation.kappa_jet(g)-baseline
                # h_kappa=-2 i/omega_g Hessian(kappa); opposite metric mode
                # conjugates h alone, leaving the complex cloud unchanged.
                dh=(-2j/omega_g*scalar_hessian(delta,omega_g,1,gamma)).conjugate()
                corrections[j,i]+=factor*np.einsum('ab,ab->',dh,raised)
        if (i+1)%8==0:print('Completed radial nodes',i+1,'/',len(radii),flush=True)
    rows=[]
    for outer,correction in zip(configurations,corrections):
        dz=kernel@correction;new=saved+dz
        rows.append(dict(analytic_kappa_outer=outer,delta_z_h=[dz.real,dz.imag],
            relative_complex_amplitude_change=float(abs(dz/saved)),
            relative_horizon_flux_change=float(abs(new/saved)**2-1),
            horizon_orbital_energy=float(data['flux']['horizon']['orbital_energy']*abs(new/saved)**2),
            maximum_absolute_projected_source_change=float(max(abs(correction)))))
    result=dict(status='isolated_kappa_boundary_flux_audit_not_full_convergence',input=path.name,
        sha256=hashlib.sha256(raw).hexdigest(),baseline_amplitude_reproduction=float(abs(recomputed/saved-1)),
        parameters=p,rows=rows,limitations=['Only dipole kappa changes; all other sectors held fixed',
        'Angular order 18 and original source cutoffs/quadrature held fixed',
        'Analytic 2000 row separates derivative implementation from boundary changes'])
    temporary=out.with_suffix('.tmp');temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(out)
    print(json.dumps(rows,indent=2),flush=True)


if __name__=='__main__':main()
