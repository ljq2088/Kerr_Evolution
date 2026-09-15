"""Independently rebuild cloud Hessian, scalar00 source and Green integral.

Shares only the physical mg=1 metric tensors and finite quadrature nodes with
production. It does NOT independently reconstruct the metric or establish
paper agreement. No fitted amplitude, ODE cloud derivative, production
Christoffel routine, scalar projection, or Green integrator is used.
"""
import hashlib,json
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.special import lpmv,gammaln,eval_legendre
from paper_leaver_cloud import LeaverThresholdCloud
from paper_leaver_radial import LeaverIngoing
from report_leaver_green_kernel import decaying_logder
from environment_radial import RadialGreen


def connection_analytic(r,t,a):
    s,c=np.sin(t),np.cos(t);sig=r*r+a*a*c*c;de=r*r-2*r+a*a
    sr,st=2*r,-2*a*a*c*s
    inv=np.zeros((4,4));inv[0,0]=-((r*r+a*a)**2-a*a*de*s*s)/(sig*de)
    inv[0,3]=inv[3,0]=-2*a*r/(sig*de);inv[3,3]=(de-a*a*s*s)/(sig*de*s*s)
    inv[1,1]=de/sig;inv[2,2]=1/sig
    dg=np.zeros((4,4,4))
    dg[1,0,0]=2/sig-2*r*sr/sig**2;dg[2,0,0]=-2*r*st/sig**2
    dg[1,0,3]=dg[1,3,0]=-2*a*s*s/sig+2*a*r*s*s*sr/sig**2
    dg[2,0,3]=dg[2,3,0]=-4*a*r*s*c/sig+2*a*r*s*s*st/sig**2
    dg[1,1,1]=(sr*de-sig*(2*r-2))/de**2;dg[2,1,1]=st/de
    dg[1,2,2]=sr;dg[2,2,2]=st
    dg[1,3,3]=2*r*s*s+2*a*a*s**4/sig-2*a*a*r*s**4*sr/sig**2
    dg[2,3,3]=2*(r*r+a*a)*s*c+8*a*a*r*s**3*c/sig-2*a*a*r*s**4*st/sig**2
    gamma=np.zeros((4,4,4))
    for k in range(4):
        for i in range(4):
            for j in range(4):
                gamma[k,i,j]=sum(inv[k,l]*(dg[i,l,j]+dg[j,l,i]-dg[l,i,j])/2 for l in range(4))
    return inv,gamma


def radial_second(cloud,r,amplitude):
    with mp.workdps(cloud.dps):
        rr=mp.mpf(str(r));d=cloud.rp-cloud.rm;x=(rr-cloud.rp)/(rr-cloud.rm)
        co=cloud.coefficients
        f=mp.polyval(list(reversed(co)),x)
        fx=mp.polyval([n*co[n] for n in range(len(co)-1,0,-1)],x)
        fxx=mp.polyval([n*(n-1)*co[n] for n in range(len(co)-1,1,-1)],x)
        xp=d/(rr-cloud.rm)**2;xpp=-2*d/(rr-cloud.rm)**3
        pref=mp.exp(-cloud.k*rr)*(rr-cloud.rm)**cloud.beta
        l=cloud.beta/(rr-cloud.rm)-cloud.k;lp=-cloud.beta/(rr-cloud.rm)**2
        return np.array([float(pref*f),float(pref*(l*f+fx*xp)),
            float(pref*((l*l+lp)*f+(2*l*xp+xpp)*fx+fxx*xp*xp))])*amplitude


def angular_second(cloud,theta):
    x=np.cos(theta);st=np.sin(theta);values=np.zeros((3,len(theta)))
    for i,b in enumerate(cloud.angular_coefficients):
        ell=1+2*i;norm=float(b)*np.sqrt((2*ell+1)/(4*np.pi)*np.exp(gammaln(ell)-gammaln(ell+2)))
        S=norm*lpmv(1,ell,x);lower=norm*lpmv(1,ell-1,x) if ell>1 else np.zeros_like(x)
        Sp=(ell*x*S-(ell+1)*lower)/st
        Spp=-x/st*Sp-(ell*(ell+1)-1/st**2)*S
        values+=np.array([S,Sp,Spp])
    return values


def scalar00_angular(a,mu,w,x):
    with mp.workdps(60):
        aa,mm,ww=map(lambda q:mp.mpf(str(q)),(a,mu,w));matrix=mp.matrix(12)
        C=lambda l:mp.mpf(l)/mp.sqrt(4*l*l-1) if l else mp.mpf(0)
        c2=aa*aa*(ww*ww-mm*mm)
        for i in range(12):
            l=2*i;matrix[i,i]=l*(l+1)-c2*(C(l+1)**2+C(l)**2)
            if i<11:matrix[i,i+1]=matrix[i+1,i]=-c2*C(l+1)*C(l+2)
        _,vectors=mp.eigsy(matrix)
        return sum(float(vectors[i,0]*mp.sign(vectors[0,0]))*np.sqrt((4*i+1)/(4*np.pi))*eval_legendre(2*i,x) for i in range(12))


def independent_source(r,theta,wt,projection,h,cloud,amplitude,angular,a):
    R,Rp,Rpp=radial_second(cloud,r,amplitude);w=float(cloud.omega);pieces=[];trace_defects=[]
    for j,t in enumerate(theta):
        S,Sp,Spp=angular[:,j];f=R*S;grad=np.array([-1j*w*f,Rp*S,R*Sp,1j*f])
        partial=np.array([[-w*w*f,-1j*w*Rp*S,-1j*w*R*Sp,w*f],
            [-1j*w*Rp*S,Rpp*S,Rp*Sp,1j*Rp*S],
            [-1j*w*R*Sp,Rp*Sp,R*Spp,1j*R*Sp],
            [w*f,1j*Rp*S,1j*R*Sp,-f]])
        inv,gamma=connection_analytic(r,t,a)
        H=partial-np.einsum('kij,k->ij',gamma,grad)
        trace_defects.append(abs(np.einsum('ij,ij->',inv,H)-float(cloud.mu)**2*f)/max(abs(float(cloud.mu)**2*f),1e-300))
        pieces.append(np.einsum('ij,ij->',inv@np.conj(h[j])@inv,H))
    sig=r*r+a*a*np.cos(theta)**2
    return 2*np.pi*np.dot(wt,projection*sig*np.array(pieces)),max(trace_defects)


def main():
    root=Path(__file__).resolve().parents[1];folder=root/'docs/environment_reproduction'
    path=folder/'fresh_20260915_L18_mg-1_sl0.json';fresh=json.loads(path.read_text());p=fresh['parameters']
    cache=json.loads((folder/'fresh_environment_audit_L18_mg1.json').read_text())['rows'][0]['metric_cache']
    metricfolder=root/'outputs/metric_cache'/cache['cache_key']
    a=p['metric']['a'];mu=p['alpha'];w=p['omega'];radii=np.array([s['r'] for s in fresh['samples']])
    weights=np.array([s['weight'] for s in fresh['samples']]);oldJ=np.array([complex(*s['source']) for s in fresh['samples']])
    theta_x,wt=np.polynomial.legendre.leggauss(p['angular_order']);theta=np.arccos(theta_x)
    cloud=LeaverThresholdCloud(mu,terms=400);energy,_=cloud.integrals(192);amplitude=1/np.sqrt(energy)
    angular=angular_second(cloud,theta);projection=scalar00_angular(a,mu,w,theta_x)
    newJ=[];source_rows=[];cache_hashes={}
    for r,old in zip(radii,oldJ):
        npz=metricfolder/(hashlib.sha256(float(r).hex().encode()).hexdigest()[:24]+'.npz')
        with np.load(npz,allow_pickle=False) as data:
            meta=json.loads(str(data['metadata']))
            assert float(data['r'])==r and meta['source_hash']==cache['source_hash']
            np.testing.assert_array_equal(meta['theta'],theta)
            assert meta['metric']['a']==a and meta['metric']['m_g']==1
            current,defect=independent_source(r,theta,wt,projection,data['h'],cloud,amplitude,angular,a)
        newJ.append(current);cache_hashes[npz.name]=hashlib.sha256(npz.read_bytes()).hexdigest()
        source_rows.append(dict(r=float(r),fresh_source=[old.real,old.imag],independent_source=[current.real,current.imag],
            relative_source_error=float(abs(current-old)/max(abs(old),1e-300)),relative_cloud_KG_trace=float(defect)))
    newJ=np.array(newJ);print('source max absolute / max old',np.max(abs(newJ-oldJ))/np.max(abs(oldJ)),flush=True)
    green=RadialGreen(a,mu,w,0,0,rmax=p['green_outer_radius'],offset=p['green_horizon_offset'],rtol=1e-11)
    series=LeaverIngoing(a,mu,w,0,0,terms=120000,dps=90);kernel=[];kernel_rows=[]
    with mp.workdps(90):
        aa,mm,ww=[mp.mpf(str(v)) for v in (a,mu,w)];match=mp.mpf(250)
        ri,rip=series.state_mp(match)
        logder,last=decaying_logder(match,aa,mm,ww,series.angular+aa*aa*ww*ww,40)
        amps=mp.lu_solve(mp.matrix([[ri,mp.conj(ri)],[rip,mp.conj(rip)]]),mp.matrix([1,logder]))
        W=(match-series.rp)*(match-series.rm)*(ri*logder-rip)
        for i,r in enumerate(radii):
            ri,_=series.state_mp(r);up=amps[0]*ri+amps[1]*mp.conj(ri);K=complex(up/W)
            kernel.append(K);reference=green.upsol.sol(r)[0]/green.w0
            kernel_rows.append(dict(r=float(r),relative_kernel_error=float(abs(K/reference-1))))
            if (i+1)%16==0:print('independent kernel',i+1,'/',len(radii),flush=True)
    kernel=np.array(kernel);oldZH=complex(*fresh['z_h']);newZH=np.sum(weights*kernel*newJ)
    kernel_only=np.sum(weights*kernel*oldJ);source_only=np.sum(weights*green.upsol.sol(radii)[0]/green.w0*newJ)
    encode=lambda z:[float(z.real),float(z.imag)]
    result=dict(status='independent_cloud_source_Green_closure_on_shared_fresh_metric',
        fresh_response_file=path.name,fresh_response_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        metric_cache=cache,metric_array_sha256=cache_hashes,source=source_rows,kernel=kernel_rows,
        source_norm_relative=float(np.max(abs(newJ-oldJ))/np.max(abs(oldJ))),
        max_kernel_relative_error=max(s['relative_kernel_error'] for s in kernel_rows),
        max_relative_cloud_KG_trace=max(s['relative_cloud_KG_trace'] for s in source_rows),
        old_z_h=encode(oldZH),independent_z_h=encode(newZH),
        complex_amplitude_relative_error=float(abs(newZH/oldZH-1)),
        horizon_flux_relative_change=float(abs(newZH/oldZH)**2-1),
        kernel_only_relative_amplitude_error=float(abs(kernel_only/oldZH-1)),
        source_only_relative_amplitude_error=float(abs(source_only/oldZH-1)),
        limitations=['Same mg=1 physical metric; no independent 10-component metric reconstruction',
          'Same finite 88-node radial quadrature, 18-node angular quadrature and source cutoffs',
          'Independent 400-term threshold cloud, unit mass integral, analytic metric derivatives, spectral field derivatives',
          'Independent 120000-term 90-digit In series and Up matching at r=250, order=40',
          'Only scalar00, not a fresh all-mode paper total; not author intermediate data'])
    if result['max_kernel_relative_error']>1e-8:
        raise RuntimeError('Independent series is not resolved on the full source grid; increase terms')
    (folder/'independent_scalar00_response_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source','kernel','metric_array_sha256','limitations']}),flush=True)


if __name__=='__main__':main()
