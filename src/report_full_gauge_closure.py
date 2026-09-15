"""Propagate every retained independently matched gauge jump to scalar00.

The change is a homogeneous pure gauge in each radial vacuum region. The
piecewise Green boundary identity therefore gives its forced scalar amplitude
without cancellation-prone numerical integration of a tiny delta-source.
"""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from report_paper_dipole_closure import old_spin1_jumps
from report_paper_kappa_tables import calculate
from paper_jump_basis import gauge_basis
from lorenz_metric import _homogeneous_radial_data
from lorenz_ghp import KerrGHP
from lorenz_mode_jet import separated_jet
from lorenz_jet import Jet
from environment_angular_diagnostic import install_dense_angular_diagnostic
from environment_source import ThresholdCloud,angular_mode
from environment_radial import RadialGreen


def enc(x):
    x=np.asarray(x,complex);return np.stack([x.real,x.imag],axis=-1).tolist()


def decode(x):
    x=np.asarray(x);return x[...,0]+1j*x[...,1]


def radial_amplitudes(r0,a,w,m,ell,spin,jumps):
    _,ru,du=_homogeneous_radial_data(spin,ell,m,a,w,r0,'Up')
    _,ri,di=_homogeneous_radial_data(spin,ell,m,a,w,r0,'In')
    A=np.array([[ru,-ri],[du,-di]]);scale=np.linalg.norm(A,axis=0)
    return np.linalg.solve(A/scale,jumps)/scale


def corrections(cloud,r0,jump_differences,inner,outer,nt=24,offset=5e-4):
    m=1;w=1/(r0**1.5+cloud.a);L=len(jump_differences)
    amps={(ell,spin):radial_amplitudes(r0,cloud.a,w,m,ell,spin,jump_differences[ell-1,k:k+2])
          for ell in range(1,L+1) for spin,k in [(-1,0),(0,2)]}
    green=RadialGreen(cloud.a,cloud.mu,cloud.omega-w,0,0,rmax=1000.,offset=1e-4,rtol=1e-11)
    xx,ww=np.polynomial.legendre.leggauss(nt);tt=np.arccos(xx)
    projection=angular_mode(tt,0,0,cloud.a**2*((cloud.omega-w)**2-cloud.mu**2))[0]
    def state(r):
        bc,index=('In',1) if r<r0 else ('Up',0);data={}
        for ell in range(1,L+1):
            for spin in (-1,0):
                _,R,Rp=_homogeneous_radial_data(spin,ell,m,cloud.a,w,r,bc)
                data[ell,spin]=np.array([R,Rp])*amps[ell,spin][index]
        Rc,Rcp=cloud.radial(r)[:,0];values=[]
        for t in tt:
            g=KerrGHP(r,t,cloud.a,omega=w,m=m,order=6)
            gc=KerrGHP(r,t,cloud.a,omega=cloud.omega,m=1,order=6)
            S,Sp,A=angular_mode(t,1,1,cloud.c2);lam=A+cloud.a**2*cloud.omega**2-2*cloud.a*cloud.omega
            field=separated_jet(gc,0,lam,Rc,Rcp,S,Sp,mass_squared=cloud.mu**2);perell=[]
            for ell in range(1,L+1):
                xi=[0 for _ in range(4)]
                for kind,spin in [('spin1',-1),('kappa',0)]:
                    for datum in (0,1):
                        basis=gauge_basis(g,ell,kind,datum,True)
                        xi=[xi[i]+data[ell,spin][datum]*basis[i] for i in range(4)]
                f=-sum(xi[i].conjugate()*gc.partial(field,i) for i in range(4))
                perell.append([f.value,f.derivative(0).value,f.derivative(0).derivative(0).value])
            values.append(perell)
        return 2*np.pi*np.einsum('k,kef->ef',ww*projection,np.array(values))
    def B(r,f):
        u,du=green.upsol.sol(r)
        return (r*r-2*r+cloud.a**2)*(u*f[:,1]-du*f[:,0])/green.w0
    end=B(outer,state(outer))-B(inner,state(inner));sides=[]
    for sign in (-1,1):
        f=state(r0+sign*offset);dr=-sign*offset
        sides.append(np.column_stack([f[:,0]+dr*f[:,1]+dr*dr*f[:,2]/2,f[:,1]+dr*f[:,2]]))
    return end+B(r0,sides[0])-B(r0,sides[1])


def main():
    p=argparse.ArgumentParser();p.add_argument('--ellmax',type=int,required=True);p.add_argument('--quadrature',type=int,default=24)
    p.add_argument('--precise-input',action='store_true')
    args=p.parse_args();L=args.ellmax
    install_dense_angular_diagnostic();Jet.coefficient_dtype=np.complex128
    root=Path(__file__).resolve().parents[1];folder=root/'docs/environment_reproduction';cloud=ThresholdCloud();r0=20.
    suffix='_precise' if args.precise_input else ''
    inputpath=folder/f'all_component_matching_L{L}_q{args.quadrature}{suffix}.json';new=json.loads(inputpath.read_text())
    for name,digest in new['implementation_sha256'].items():
        if hashlib.sha256((root/'src'/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('Matching input implementation changed: '+name)
    from source_provenance import local_dependency_hashes,samples_hash
    provenance=dict(files=local_dependency_hashes(root/'src',['report_paper_dipole_closure','report_paper_kappa_tables']),angular='DenseRealHarmonic',jet_dtype='complex128',r0=r0,a=cloud.a,m=1)
    digest=hashlib.sha256(json.dumps(provenance,sort_keys=True).encode()).hexdigest()
    oldcache=folder/f'full_gauge_legacy_jump_inputs_{digest[:12]}.json'
    olddata=json.loads(oldcache.read_text()) if oldcache.exists() else {'provenance':provenance,'rows':[]}
    if olddata['provenance']!=provenance:raise ValueError('Legacy basis inputs changed')
    if olddata['rows'] and olddata.get('rows_sha256')!=samples_hash(olddata['rows']):raise ValueError('Legacy basis input checksum changed')
    for ell in range(len(olddata['rows'])+1,L+1):
        spin,error=old_spin1_jumps(r0,cloud.a,1,ell)
        k0,k1,_=calculate(cloud.a,r0,1,ell)
        olddata['rows'].append(dict(ell=ell,jumps=enc(np.r_[spin,k0,k1]),basis_fit_relative_error=error))
        olddata['rows_sha256']=samples_hash(olddata['rows'])
        oldcache.write_text(json.dumps(olddata,indent=2)+'\n');print('legacy jumps',ell,error,flush=True)
    Jet.coefficient_dtype=np.clongdouble
    old=np.array([decode(row['jumps']) for row in olddata['rows'][:L]])
    current=decode(new['solution']).reshape(L,4);difference=current-old
    fresh=folder/'fresh_20260915_L18_mg-1_sl0.json';baseline=json.loads(fresh.read_text());Z=complex(*baseline['z_h'])
    pars=baseline['parameters'];dz=corrections(cloud,r0,difference,pars['source_panels'][0],pars['source_panels'][-1])
    total=np.sum(dz);rows=[]
    for ell in range(1,L+1):
        rows.append(dict(ell=ell,old_jumps=enc(old[ell-1]),new_jumps=enc(current[ell-1]),jump_difference=enc(difference[ell-1]),
            delta_z_h=enc(dz[ell-1]),cumulative_flux_relative_change=float(abs(1+np.sum(dz[:ell])/Z)**2-1)))
    d=dict(status='all_retained_gauge_jump_corrections_propagated_not_full_paper_reproduction',L=L,
        input_matching_sha256=hashlib.sha256(inputpath.read_bytes()).hexdigest(),
        matching_input=str(inputpath.relative_to(root)),legacy_inputs=str(oldcache.relative_to(root)),
        legacy_inputs_sha256=hashlib.sha256(oldcache.read_bytes()).hexdigest(),
        implementation_sha256=local_dependency_hashes(root/'src',['report_full_gauge_closure']),
        input_response_sha256=hashlib.sha256(fresh.read_bytes()).hexdigest(),baseline_z_h=enc(Z),delta_z_h=enc(total),
        relative_flux_change=float(abs(1+total/Z)**2-1),rows=rows,
        limitations=['Only the finite L gauge-amplitude difference is changed; curvature and trace inputs are shared',
          'Baseline has L18 while only the first L sectors receive independently matched gauge corrections',
          'Held-out tensor matching and L convergence must be inspected separately',
          'No author metric arrays and no fitted environmental flux'])
    out=folder/f'full_gauge_closure_L{L}{suffix}.json';out.write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(d),flush=True)


if __name__=='__main__':main()
