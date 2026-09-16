"""Fresh dipole horizon-source endpoint test, not fixed-source Green precision.

A pure gauge dipole scalar obeys L f=J separately on each vacuum side. For
lowering the inner source cutoff only, delta Z_H=B(old)-B(new), with
B=Delta*(R_up*fprime-R_upprime*f)/W. Orbit terms cancel in this difference.
Trace/kappa diagnostic integration is extended below all new evaluation radii;
no RadialGreen dense-output extrapolation below its horizon start is allowed.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
import lorenz_kappa
from lorenz_metric import spin0_metric,spin1_metric
from lorenz_jet import Jet
from lorenz_ghp import KerrGHP
from lorenz_mode_jet import separated_jet
from environment_source import ThresholdCloud,angular_mode
from environment_radial import RadialGreen
from environment_cloud import mode_flux
from source_provenance import source_fingerprint,local_dependency_hashes
ROOT=Path(__file__).resolve().parents[1]
def enc(z):
    z=np.asarray(z,complex);return np.stack([z.real,z.imag],axis=-1).tolist()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    pa=argparse.ArgumentParser();pa.add_argument("--baseline",required=True);pa.add_argument('--output',required=True);pa.add_argument('--angular-order',type=int,default=18);pa.add_argument('--offset',type=float,default=1e-7);pa.add_argument('--rtol',type=float,default=1e-12);pa.add_argument('--extended',action='store_true');args=pa.parse_args()
    inp=ROOT/args.baseline;old=json.loads(inp.read_text());par=old['parameters'];cloud=ThresholdCloud(alpha=par['alpha']);r0=par['metric']['orbital_radius'];omega=par['omega'];scalar_ell=par['scalar_ell'];scalar_m=par['scalar_m'];metric_m=par['metric']['m_g']
    if abs(metric_m)!=1 or scalar_m!=metric_m+1 or par['cloud_mass']!=1:
        raise ValueError('Requires a unit-mass threshold cloud and metric dipole mg=+/-1')
    if abs(par['source_panels'][0]-cloud.rp-5e-4)>1e-12:
        raise ValueError('Baseline source cutoff must be rplus+5e-4')
    out=Path(args.output)
    if out.exists():raise FileExistsError(out)
    if args.extended:Jet.coefficient_dtype=np.clongdouble
    green=RadialGreen(cloud.a,cloud.mu,omega,scalar_ell,scalar_m,rmax=par['green_outer_radius'],offset=args.offset,rtol=args.rtol,infinity_method=par.get('infinity_method','series'))
    rnodes=np.array([s['r'] for s in old['samples']]);wgt=np.array([s['weight'] for s in old['samples']]);J=np.array([complex(*s['source']) for s in old['samples']]);Z=np.dot(wgt,green.upsol.sol(rnodes)[0]*J)/green.w0
    Zbase=complex(*old['z_h']);base_reint=float(abs(Z-Zbase)/abs(Zbase))
    if base_reint>1e-8:raise ValueError('Extended scalar Green fails baseline amplitude check')
    xx,ww=np.polynomial.legendre.leggauss(args.angular_order);theta=np.arccos(xx);Sout=angular_mode(theta,scalar_ell,scalar_m,cloud.a**2*(omega**2-cloud.mu**2))[0]
    def state(r):
        if r<max(cloud.rmin,green.rmin):raise ValueError('Evaluation outside solved cloud/Green domain')
        R,Rp=cloud.radial(r)[:,0];values=[]
        for t in theta:
            g,one=spin1_metric(r,t,r0,cloud.a,1,1,6,return_vector=True,full_current=True)
            _,zero=spin0_metric(r,t,r0,cloud.a,1,1,6,return_vector=True)
            gc=KerrGHP(r,t,cloud.a,omega=cloud.omega,m=1,order=6)
            ss,sp,A=angular_mode(t,1,1,cloud.c2);lam=A+cloud.a**2*cloud.omega**2-2*cloud.a*cloud.omega
            field=separated_jet(gc,0,lam,R,Rp,ss,sp,mass_squared=cloud.mu**2)
            fv=[]
            for vec in [zero,one]:
                lie=[((-1j*v/g.omega).conjugate() if metric_m<0 else -1j*v/g.omega) for v in vec]
                f=sum(g.inv[i][j]*lie[j]*gc.partial(field,i) for i in range(4) for j in range(4))
                fv.append([f.value,f.derivative(0).value])
            values.append(fv)
        return 2*np.pi*np.einsum('k,ksd->sd',ww*Sout,np.array(values))
    def boundary(r,st):
        u,up=green.upsol.sol(r);return (r*r-2*r+cloud.a**2)*(u*st[:,1]-up*st[:,0])/green.w0
    # Record old-domain value at its safe cutoff before modifying only local bindings.
    initial_offset=5e-4;initial_r=cloud.rp+initial_offset;oldstate=state(initial_r);oldB=boundary(initial_r,oldstate)
    original_class=lorenz_kappa.RadialGreen
    def extended_radial(*a,**kw):
        kw['offset']=args.offset;kw['rtol']=args.rtol;return original_class(*a,**kw)
    lorenz_kappa.RadialGreen=extended_radial;lorenz_kappa.scalar_resolvent.cache_clear()
    eps=[5e-4,2.5e-4,5e-5,2.5e-5,5e-6,2.5e-6,1.2e-6]
    result=dict(status='sampling_horizon_dipole_cutoff',baseline=str(inp.relative_to(ROOT)),baseline_sha256=sha(inp),baseline_z_h=enc(Zbase),baseline_flux=old['flux']['horizon']['orbital_energy'],baseline_source_provenance=old.get('source_provenance'),current_source_fingerprint=source_fingerprint(),parameters=dict(a=cloud.a,rp=r0,alpha=cloud.mu,omega=omega,scalar_ell=scalar_ell,scalar_m=scalar_m,metric_m=metric_m,cloud_rmin=cloud.rmin,scalar_green_rmin=green.rmin,trace_kappa_rmin=cloud.rp+args.offset,offset=args.offset,rtol=args.rtol,angular_order=args.angular_order,jet_dtype=str(Jet.coefficient_dtype)),baseline_extended_green_reintegration_relative=base_reint,original_safe_cutoff_state=enc(oldstate),original_safe_cutoff_boundary_by_sector=enc(oldB),formula='delta Z_H(old->new)=sum_s[B_s(old)-B_s(new)], B_s=Delta*(Rup*f_sprime-Rupprime*f_s)/W. Orbit terms cancel because only inner cutoff moves.',implementation_sha256=local_dependency_hashes(ROOT/'src',['report_field_horizon_cutoff']),rows=[],limitations=['Stored baseline is historical and not relabelled as freshly sampled; fresh current-code dipole only; fullL18 baseline corrected by this dipole shell contribution, not an all-ell source extension.','Current cloud normalization unchanged; lowest radius remains above cloud integration start.','Gauge identity is valid within inner vacuum region; no orbit distribution is added.','Trace/kappa and massive scalar domains explicitly extended; no extrapolation below their starting radius.'])
    def save():
        tmp=out.with_suffix('.tmp');tmp.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');tmp.replace(out)
    start=time.perf_counter();Bref=None
    for cut in eps:
        rr=cloud.rp+cut;st=state(rr);B=boundary(rr,st)
        if Bref is None:Bref=B.copy()
        delta=Bref-B;Zn=Zbase+sum(delta);F=mode_flux(omega,scalar_m,cloud.omega,1,cloud.mu,cloud.a,0j,Zn)['horizon']['orbital_energy'];F0=result['baseline_flux']
        result['rows'].append(dict(offset=cut,r=rr,state_by_spin0_spin1=enc(st),boundary_by_spin0_spin1=enc(B),delta_Z_by_spin0_spin1=enc(delta),delta_Z=enc(sum(delta)),corrected_full_baseline_Z=enc(Zn),corrected_full_baseline_flux=F,relative_flux_change=F/F0-1,relative_amplitude_correction=float(abs(sum(delta))/abs(Zbase))))
        save();print('cutoff',cut,'deltaZ',sum(delta),'dF/F',F/F0-1,flush=True)
    result['old_vs_extended_safe_cutoff_B_difference']=enc(Bref-oldB)
    # B(x)=B0+c*x+d*x^(1-2i gamma)+higher powers, from the two Frobenius phases.
    gamma=2*cloud.rp*(omega-scalar_m*cloud.a/(2*cloud.rp))/(2*np.sqrt(1-cloud.a**2));offsets=np.array(eps);vals=np.array([complex(*row['boundary_by_spin0_spin1'][0])+complex(*row['boundary_by_spin0_spin1'][1]) for row in result['rows']]);fits=[]
    for degree in [1,2]:
        for count in [5,7]:
            x=offsets[-count:];y=vals[-count:];scale=max(x);t=x/scale
            cols=[np.ones(count,complex)]
            for n in range(1,degree+1):cols.extend([t**n,t**n*np.exp(-2j*gamma*np.log(x))])
            A=np.array(cols).T;co,_,rank,_=np.linalg.lstsq(A,y,rcond=None);B0=co[0];delta=sum(Bref)-B0;Zn=Zbase+delta;F=mode_flux(omega,scalar_m,cloud.omega,1,cloud.mu,cloud.a,0j,Zn)['horizon']['orbital_energy']
            fits.append(dict(degree=degree,sample_count=count,rank=int(rank),condition=float(np.linalg.cond(A)),residual_max=float(np.max(abs(A@co-y))),horizon_B=enc(B0),estimated_total_omitted_Z=enc(delta),relative_flux_change=F/result['baseline_flux']-1))
    result.update(status='completed_fresh_dipole_horizon_cutoff',frobenius_gamma=gamma,fits=fits,accepted_overdetermined_fits=[f for f in fits if f['sample_count']>1+2*f['degree']],elapsed_seconds=time.perf_counter()-start)
    save()
if __name__=='__main__':main()
