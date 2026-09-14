"""Independent boundary identity for the actual scalar00 dipole source.

In each vacuum region the metric dipole is pure gauge. Its induced scalar f
obeys L f=J. Green's identity converts the source integral to endpoint and
one-sided orbit terms; the latter MUST remain when f is discontinuous.
"""
import hashlib,json
from pathlib import Path
import numpy as np
from environment_source import ThresholdCloud,angular_mode
from environment_radial import RadialGreen
from lorenz_ghp import KerrGHP
from lorenz_mode_jet import separated_jet
from lorenz_metric import spin1_metric,spin0_metric


def main():
    folder=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    path=folder/'forced_mode_nr16_nt18_L1_mg-1_sl0_inner0.0005_outer320_log_h64.json'
    raw=path.read_bytes();d=json.loads(raw);p=d['parameters']
    assert p['metric']['m_g']==-1 and p['metric']['ellmax']==1 and p['scalar_ell']==p['scalar_m']==0
    cloud=ThresholdCloud(alpha=p['alpha']);r0=p['metric']['orbital_radius']
    green=RadialGreen(cloud.a,cloud.mu,p['omega'],0,0,rmax=p['green_outer_radius'],
                      offset=p['green_horizon_offset'],rtol=1e-11)
    x,w=np.polynomial.legendre.leggauss(p['angular_order']);theta=np.arccos(x)
    angular=angular_mode(theta,0,0,cloud.a**2*(p['omega']**2-cloud.mu**2))[0]
    samples=[]
    def state(r):
        R,Rp=cloud.radial(r)[:,0];values=[]
        for t in theta:
            g,one=spin1_metric(r,t,r0,cloud.a,1,1,6,return_vector=True,full_current=True)
            _,zero=spin0_metric(r,t,r0,cloud.a,1,1,6,return_vector=True)
            minus_lie_vector=[(-1j*(one[i]+zero[i])/g.omega).conjugate() for i in range(4)]
            gc=KerrGHP(r,t,cloud.a,omega=cloud.omega,m=1,order=6)
            S,Sp,A=angular_mode(t,1,1,cloud.c2)
            eigenvalue=A+cloud.a**2*cloud.omega**2-2*cloud.a*cloud.omega
            field=separated_jet(gc,0,eigenvalue,R,Rp,S,Sp,mass_squared=cloud.mu**2)
            f=sum(g.inv[i][j]*minus_lie_vector[j]*gc.partial(field,i) for i in range(4) for j in range(4))
            values.append([f.value,f.derivative(0).value,f.derivative(0).derivative(0).value])
        projected=2*np.pi*(w*angular)@np.array(values)
        samples.append(dict(r=r,field_and_derivatives=enc(projected)))
        print('Gauge scalar sampled',r,flush=True)
        return projected
    def boundary(r,st):
        up,dup=green.upsol.sol(r)
        return (r*r-2*r+cloud.a**2)*(up*st[1]-dup*st[0])/green.w0
    inner,outer=p['source_panels'][0],p['source_panels'][-1]
    endpoint=boundary(outer,state(outer))-boundary(inner,state(inner))
    saved=complex(*d['z_h']);rows=[]
    for epsilon in (.001,.0005):
        limits=[]
        for sign in (-1,1):
            st=state(r0+sign*epsilon);dr=-sign*epsilon
            limits.append(np.array([st[0]+dr*st[1]+dr*dr*st[2]/2,st[1]+dr*st[2]]))
        orbit=boundary(r0,limits[0])-boundary(r0,limits[1])
        predicted=endpoint+orbit
        rows.append(dict(epsilon=epsilon,endpoint_term=enc(endpoint),orbit_jump_term=enc(orbit),
            predicted_z_h=enc(predicted),relative_amplitude_error=float(abs(predicted/saved-1)),
            relative_flux_change=float(abs(predicted/saved)**2-1)))
    result=dict(status='dipole_green_identity_not_full_metric_validation',input=path.name,
        sha256=hashlib.sha256(raw).hexdigest(),saved_z_h=enc(saved),rows=rows,samples=samples,
        limitations=['Pure gauge identity applies only within each vacuum region',
                     'Orbit terms are not extra physical delta sources or an authorized correction',
                     'Same source cutoffs and angular quadrature; no full metric or gauge equivalence proof'])
    out=folder/'dipole_scalar00_green_identity.json';out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(rows,indent=2),flush=True)


def enc(a):
    a=np.asarray(a);return np.stack([a.real,a.imag],axis=-1).tolist()


if __name__=='__main__':main()
