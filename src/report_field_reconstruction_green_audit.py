"""Independent fixed-Gauss Green reconstruction and angular/Fourier field audit."""
import hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.interpolate import BarycentricInterpolator
from scipy.special import obl_ang1,pro_ang1,sph_harm
from environment_wake import EnvironmentalWake
from environment_source import angular_mode
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'docs/environment_reproduction'
def enc(v):
    a=np.asarray(v,complex);return np.stack([a.real,a.imag],axis=-1).tolist()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def independent_angular(theta,ell,m,c2):
    mm=abs(m);f=obl_ang1 if c2>=0 else pro_ang1;c=np.sqrt(abs(c2));x,w=np.polynomial.legendre.leggauss(96)
    raw=f(mm,ell,c,x)[0];base=sph_harm(mm,ell,np.zeros(len(x)),np.arccos(x)).real
    sign=np.sign(np.dot(w,raw*base));scale=sign/np.sqrt(2*np.pi*np.dot(w,raw*raw))
    return scale*f(mm,ell,c,np.cos(theta))[0]*((-1)**mm if m<0 else 1)
def direct_field(g,data,radii,n):
    """Integrate independently constructed barycentric source polynomials.

    Fixed Gauss quadrature splits every panel at each requested radius. No
    SampledResponse prefix/suffix or its ODE dense output is used.
    """
    par=data['parameters'];nodes=np.array([r['r'] for r in data['samples']]);J=np.array([complex(*r['source']) for r in data['samples']]);panels=par['source_panels'];q,w=np.polynomial.legendre.leggauss(n);A=np.zeros(len(radii),complex);B=A.copy()
    for panel,(lo,hi) in enumerate(zip(panels[:-1],panels[1:])):
        mask=(nodes>lo)&(nodes<hi);log=(panel==0 and par.get('horizon_log_first_panel',False));zl,zh=np.log([lo-g.rp,hi-g.rp]) if log else (lo,hi)
        xn,wn=np.polynomial.legendre.leggauss(sum(mask));bary=(-1.)**np.arange(len(xn))*np.sqrt((1-xn*xn)*wn);poly=BarycentricInterpolator(xn,J[mask],wi=bary)
        def integral(a,b,branch):
            if a>=b:return 0j
            za,zb=np.log([a-g.rp,b-g.rp]) if log else (a,b);z=(za+zb)/2+(zb-za)*q/2;r=g.rp+np.exp(z) if log else z;jac=np.exp(z) if log else 1.
            source=poly(2*(z-zl)/(zh-zl)-1);sol=g.insol if branch=='In' else g.upsol
            return (zb-za)/2*np.dot(w,jac*sol.sol(r)[0]*source/g.w0)
        for i,r in enumerate(radii):
            if r>lo:A[i]+=integral(lo,min(r,hi),'In')
            if r<hi:B[i]+=integral(max(r,lo),hi,'Up')
    u,du=g.insol.sol(radii);v,dv=g.upsol.sol(radii)
    return np.array([v*A+u*B,dv*A+du*B])
def main():
    out=D/'field_reconstruction_green_audit_20260916.json';start=time.perf_counter();manifest={};rows=[];groups=[]
    def load(p):manifest[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
    fig=load(D/'figure5_threshold_wakes.json');cov20=load(D/'flux_coverage_L18_nt18_i12_h12_f12_m0go4000_m1go4000_m2go4000.json')
    sources20={(r['ell'],r['m']):D/r['file'] for r in cov20['field']['modes'] if 'flux' in r}
    sources416={}
    for item in fig['datasets'][0]['inputs']:
        p=ROOT/item['file'];d=json.loads(p.read_text());sources416[(d['parameters']['scalar_ell'],d['parameters']['scalar_m'])]=p
    plan=[(20.,sources20,[(2,-2),(2,0),(2,2),(3,1),(3,3)]),(41.6,sources416,[(2,-2),(2,0),(2,2),(3,1),(3,3),(4,2),(6,2)])]
    result=dict(status='in_progress',method='Independent barycentric source interpolation and fixed split Gauss64/128 Green integral; independent scipy spheroidal value/phase; saved full-field Fourier decomposition',rows=rows,groups=groups)
    def save():result['inputs_sha256']=manifest;out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for r0,sources,keys in plan:
        data=[load(sources[key]) for key in keys];wake=EnvironmentalWake(data,ellmax=6,allow_partial=True,allow_mixed_discretization=True);rp=1+np.sqrt(1-wake.parameters['a']**2);radii=np.array([rp+.00025,1.5430278953672374,3.,r0,100.,200.]);reference_fields={}
        for d in data:
            p=d['parameters'];key=(p['scalar_ell'],p['scalar_m']);omega,response=wake.modes[key];actual=response.evaluate(radii);v64=direct_field(response.green,d,radii,64);v128=direct_field(response.green,d,radii,128);scale=np.maximum(np.max(abs(v128),axis=1)[:,None]*1e-12,abs(v128))
            c2=wake.parameters['a']**2*(omega**2-wake.parameters['alpha']**2);theta=np.array([.4,1.1,np.pi/2,2.2]);sv=angular_mode(theta,*key,c2)[0];ind=independent_angular(theta,*key,c2)
            saved_z=complex(*d['z_h']);orb=next(v for v in d['radial_response_at_panel_boundaries'] if v['r']==r0)
            row=dict(orbit=r0,ell=key[0],m=key[1],input=str(sources[key].relative_to(ROOT)),radii=radii.tolist(),sampled_response=enc(actual),direct_gauss64=enc(v64),direct_gauss128=enc(v128),relative64_to128=(abs(v64-v128)/scale).tolist(),relative_sampled_to128=(abs(actual-v128)/scale).tolist(),saved_gauss_at_orbit=orb,relative_particle_field_saved_to_direct=abs(complex(*orb['field'])-v128[0,3])/max(abs(v128[0,3]),1e-30),angular_max_absolute_error=float(np.max(abs(sv-ind))),horizon_amplitude_from_flux=np.sqrt(abs(d['flux']['horizon']['orbital_energy']/(4*rp*(omega-(wake.omega_c))*(omega-key[1]*response.green.oh)))) if omega!=wake.omega_c else None,saved_ZH_absolute=abs(saved_z),continuous_ZH=enc(response.horizon_coefficient),relative_continuous_ZH_to_saved=abs(response.horizon_coefficient-saved_z)/max(abs(saved_z),1e-300),sourcefree_identity_absolute_error=float(abs(actual[0,0]-response.green.insol.sol(radii[0])[0]*response.horizon_coefficient)))
            rows.append(row);reference_fields[key]=v128[0];save();print('audited',r0,key,'response maxrel',np.max(abs(actual-v128)/scale),'angular',row['angular_max_absolute_error'],flush=True)
        # Independent angular factors and direct Green values at three spacetime points.
        point_r=radii[[1,3,4]];theta=np.array([np.pi/2,1.1,.7]);phi=np.array([.3,-.4,1.2]);event_time=17.3;manual=np.zeros(3,complex)
        for key,(omega,response) in wake.modes.items():
            c2=wake.parameters['a']**2*(omega**2-wake.parameters['alpha']**2);manual+=reference_fields[key][[1,3,4]]*independent_angular(theta,*key,c2)*np.exp(1j*(key[1]*phi-omega*event_time))
        assembled=wake.evaluate(point_r,theta,phi,time=event_time);dt=3.7;hel=wake.evaluate(point_r,theta,phi+wake.omega_p*dt,time=event_time+dt);expected=assembled*np.exp(-1j*(wake.omega_c-wake.omega_p)*dt)
        groups.append(dict(orbit=r0,modes=keys,manual_independent_field=enc(manual),environmental_wake=enc(assembled),relative_assembly_error=(abs(manual-assembled)/np.maximum(abs(manual),1e-30)).tolist(),helical_covariance_max_absolute_error=float(np.max(abs(hel-expected))),limitation='Selected partial modes validate assembly; missing modes are not set to zero in the physical Fig5 report.'))
        save()
    # Fourier analysis of already saved complete complex fields is exact on this grid.
    npz=D/'figure5_threshold_wakes.npz';manifest[str(npz.relative_to(ROOT))]=sha(npz);arr=np.load(npz);phi=(arr['angular_edges'][:-1]+arr['angular_edges'][1:])/2;fourier=[]
    for item,tag in zip(fig['datasets'],['416','418']):
        field=arr['field_'+tag];rr=np.sqrt(arr['radial_edges_'+tag][:-1]*arr['radial_edges_'+tag][1:]);idx=np.unravel_index(np.argmax(abs(field)),field.shape);f=field[idx[0]];coef={m:np.dot(f,np.exp(-1j*m*phi))/len(phi) for m in range(-12,13)};rms=np.sqrt(np.mean(abs(f)**2));recon=sum(v*np.exp(1j*m*phi) for m,v in coef.items());m2=[]
        for inp in item['inputs']:
            p=ROOT/inp['file'];d=json.loads(p.read_text());par=d['parameters']
            if par['scalar_m']==2:
                manifest[str(p.relative_to(ROOT))]=sha(p);ang=angular_mode(np.pi/2,par['scalar_ell'],2,par['metric']['a']**2*(par['omega']**2-par['alpha']**2))[0];m2.append(dict(ell=par['scalar_ell'],horizon_equatorial_amplitude_per_epsilon=abs(complex(*d['z_h'])*ang)/par['alpha']**3))
        a2=next(v['horizon_equatorial_amplitude_per_epsilon'] for v in m2 if v['ell']==2);lower=a2-sum(v['horizon_equatorial_amplitude_per_epsilon'] for v in m2 if v['ell']!=2)
        fourier.append(dict(orbit=item['parameters']['r0'],peak_r=float(rr[idx[0]]),peak_phi=float(phi[idx[1]]),maximum=float(abs(field[idx])),angular_rms=float(rms),dominant_Fourier_modes=[dict(m=m,coefficient=enc(v),absolute=float(abs(v))) for m,v in sorted(coef.items(),key=lambda p:abs(p[1]),reverse=True)[:6]],fourier_reconstruction_max_error=float(np.max(abs(f-recon))),parseval_absolute_error=float(abs(rms*rms-sum(abs(v)**2 for v in coef.values()))),m2_horizon_mode_amplitudes=m2,phase_independent_m2_horizon_lower_bound=float(lower),interpretation='Changing relative m phases preserves angular RMS, so cannot make the maximum smaller than that RMS. Even arbitrary phases between the retained ell within m2 cannot reduce its horizon coefficient below the triangle lower bound. This is conditional on the saved physical mode amplitudes and the same field definition, not a statement that an original figure colorbar is its global maximum.'))
    result.update(status='completed_local_field_reconstruction_audit',fourier_full_saved_field=fourier,elapsed_seconds=time.perf_counter()-start,implementation_sha256=sha(Path(__file__)),limitations=['Same saved source functions and homogeneous radial solutions are shared; no independent metric/source derivation.','Agreement of interpolation/integration paths does not certify source interpolation between its original sample nodes.','Selected rp20/rp41.6 modes, not a recomputed complete field.','No mode deletion, complex phase fit, or amplitude correction applied.']);save()
if __name__=='__main__':main()
