"""Propagate independently matched dipole jump corrections to scalar00 Z_H.

Uses pure-gauge Green identity for the homogeneous correction. No fitting of
any parameter to the environmental flux or paper curve is performed.
"""
import json
from pathlib import Path
import numpy as np
from paper_jump_basis import gauge_basis,radial_jet,angular_jet,operators
from paper_sourced_matching import orbit
from report_paper_kappa_tables import calculate
from lorenz_metric import spin1_metric,_homogeneous_radial_data
from lorenz_ghp import KerrGHP
from lorenz_mode_jet import separated_jet
from environment_source import ThresholdCloud,angular_mode
from environment_radial import RadialGreen


def old_spin1_jumps(r0,a,m,ell,eps=1e-4):
    rows=[];vals=[]
    for t in (.6,1.1,1.5,2.,2.6):
        g=KerrGHP(r0,t,a,omega=m/(r0**1.5+a),m=m,order=6)
        basis=[gauge_basis(g,ell,'spin1',d,True) for d in (0,1)]
        rows.extend([[v[i].value for v in basis] for i in range(4)])
        vjump=np.zeros(4,complex)
        for side in (-1,1):
            gs,raw=spin1_metric(r0+side*eps,t,r0,a,ell,m,order=8,return_vector=True,full_current=True)
            v=[1j/gs.omega*sum(gs.inv[i][k]*raw[k] for k in range(4)) for i in range(4)]
            dr=-side*eps
            vjump+=side*np.array([u.value+dr*u.derivative(0).value+dr*dr/2*u.derivative(0).derivative(0).value for u in v])
        vals.extend(vjump)
    A=np.array(rows);b=np.array(vals);x=np.linalg.lstsq(A,b,rcond=None)[0]
    return x,float(np.linalg.norm(A@x-b)/np.linalg.norm(b))


def amplitudes(r0,a,w,m,spin,jumps):
    _,ru,du=_homogeneous_radial_data(spin,1,m,a,w,r0,'Up')
    _,ri,di=_homogeneous_radial_data(spin,1,m,a,w,r0,'In')
    return np.linalg.solve(np.array([[ru,-ri],[du,-di]]),jumps)


def response(cloud,r0,delta,inner,outer,nt=18,offset=5e-4):
    m=1;w=1/(r0**1.5+cloud.a)
    amps={-1:amplitudes(r0,cloud.a,w,m,-1,delta[:2]),0:amplitudes(r0,cloud.a,w,m,0,delta[2:])}
    green=RadialGreen(cloud.a,cloud.mu,cloud.omega-w,0,0,rmax=1000.,offset=1e-4,rtol=1e-11)
    xx,ww=np.polynomial.legendre.leggauss(nt);tt=np.arccos(xx)
    angular=angular_mode(tt,0,0,cloud.a**2*((cloud.omega-w)**2-cloud.mu**2))[0]
    def state(r):
        bc,idx=('In',1) if r<r0 else ('Up',0)
        data={}
        for spin in (-1,0):
            _,R,Rp=_homogeneous_radial_data(spin,1,m,cloud.a,w,r,bc)
            data[spin]=np.array([R,Rp])*amps[spin][idx]
        Rc,Rcp=cloud.radial(r)[:,0];vals=[]
        for t in tt:
            g=KerrGHP(r,t,cloud.a,omega=w,m=m,order=6)
            xi=[0 for _ in range(4)]
            for kind,spin in (('spin1',-1),('kappa',0)):
                for datum in (0,1):
                    basis=gauge_basis(g,1,kind,datum,True)
                    xi=[xi[i]+data[spin][datum]*basis[i] for i in range(4)]
            gc=KerrGHP(r,t,cloud.a,omega=cloud.omega,m=1,order=6)
            S,Sp,A=angular_mode(t,1,1,cloud.c2)
            lam=A+cloud.a**2*cloud.omega**2-2*cloud.a*cloud.omega
            field=separated_jet(gc,0,lam,Rc,Rcp,S,Sp,mass_squared=cloud.mu**2)
            # Physical scalar mg=-1 uses conjugated metric gauge vector only.
            f=-sum(xi[i].conjugate()*gc.partial(field,i) for i in range(4))
            vals.append([f.value,f.derivative(0).value,f.derivative(0).derivative(0).value])
        return 2*np.pi*(ww*angular)@np.array(vals)
    def B(r,f):
        u,du=green.upsol.sol(r)
        return (r*r-2*r+cloud.a**2)*(u*f[1]-du*f[0])/green.w0
    end=B(outer,state(outer))-B(inner,state(inner));sides=[]
    for sign in (-1,1):
        f=state(r0+sign*offset);dr=-sign*offset
        sides.append(np.array([f[0]+dr*f[1]+dr*dr*f[2]/2,f[1]+dr*f[2]]))
    return end+B(r0,sides[0])-B(r0,sides[1])


def enc(x):
    z=np.asarray(x,complex);return np.stack([z.real,z.imag],axis=-1).tolist()


def main():
    folder=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    cloud=ThresholdCloud();r0=20.;spin,fiterror=old_spin1_jumps(r0,cloud.a,1,1)
    k0,k1,_=calculate(cloud.a,r0,1,1);old=np.r_[spin,k0,k1]
    original=json.loads((folder/'forced_mode_nr16_nt18_L1_mg-1_sl0_inner0.0005_outer320_log_h64.json').read_text())
    Z=complex(*original['z_h']);p=original['parameters'];rows=[]
    for L in (4,6):
        d=json.loads((folder/f'paper_sourced_jumps_r20_m1_L{L}_extended.json').read_text())
        v=np.array(d['solution']);new=v[:4,0]+1j*v[:4,1];delta=new-old
        dz=response(cloud,r0,delta,p['source_panels'][0],p['source_panels'][-1])
        rows.append(dict(L=L,matched_jumps=enc(new),legacy_jumps=enc(old),jump_difference=enc(delta),
             relative_jump_difference=enc(delta/old),delta_z_h=enc(dz),
             relative_flux_change=float(abs((Z+dz)/Z)**2-1)))
        print(json.dumps(rows[-1]),flush=True)
    result=dict(status='independent_local_matching_correction_propagated_to_dipole_flux',
      old_vector_fit_relative_residual=fiterror,old_z_h=enc(Z),rows=rows,
      limitations=['Only ell=1 homogeneous gauge correction propagated',
        'New matching uses extended jets with double angular/radial input and QR',
        'Not a complete replacement metric or a new total flux'])
    (folder/'paper_dipole_closure.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
