"""Isolated multidomain Chebyshev BVP test on frozen Fig.10 sources.
No production solver, caches or active run files are changed.
Solve (Delta R')'+V R=J using z=log(r-r+), piecewise Lobatto
collocation, continuity of R and R', and analytic finite-end Robin data.
This tests radial integration, not metric reconstruction/source correctness.
"""
from pathlib import Path
import json, hashlib, sys, time
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from numpy.polynomial import legendre as L, chebyshev as C
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from environment_cloud import angular_eigenvalue, radial_coefficients, mode_flux
from environment_radial import horizon_series
from report_li_order4_boundary_control import formal_boundary
from li_order4_green import LiOrder4Green
from environment_response import SampledResponse
OUT=Path(__file__).resolve().parent

def grid(n):
    t=-np.cos(np.arange(n+1)*np.pi/n)
    w=(-1.)**np.arange(n+1);w[[0,-1]]*=.5
    dx=t[:,None]-t[None,:];np.fill_diagonal(dx,1.)
    D=(w[None,:]/w[:,None])/dx;np.fill_diagonal(D,0.)
    D[np.diag_indices(n+1)]=-D.sum(axis=1)
    assert np.max(abs(D@t**3-3*t**2))<1e-9
    return t,D

def source_parts(d,rp):
    panels=np.array([rp+5e-4,3,6,12,20,40,80,160,320.])
    r=np.array(d['radii']);J=np.array([complex(*q) for q in d['source']])
    parts=[]
    for j,(a,b) in enumerate(zip(panels[:-1],panels[1:])):
        mask=(r>a)&(r<b);z=np.log(r[mask]-rp) if j==0 else r[mask]
        lo,hi=np.log([a-rp,b-rp]) if j==0 else (a,b)
        x=2*(z-lo)/(hi-lo)-1
        assert np.allclose(x,L.leggauss(len(x))[0],atol=2e-9,rtol=0)
        parts.append((lo,hi,L.legfit(x,J[mask],len(x)-1),j==0))
    return panels,r,J,parts

def source_eval(rad,index,parts,rp):
    if index is None:return np.zeros_like(rad,dtype=complex)
    lo,hi,c,log=parts[index];z=np.log(rad-rp) if log else rad
    return L.legval(2*(z-lo)/(hi-lo)-1,c)

def spectral(d,n,parts,panels,outer=1000.):
    a=.88;mu=.3;omega=d['omega'];m=d['m'];ell=d['ell']
    assert 0 < omega < mu, 'This initial trial is restricted to bound channels.'
    rp=1+np.sqrt(1-a*a);lam=angular_eigenvalue(ell,m,a*a*(omega*omega-mu*mu))
    edges=np.r_[rp+1e-4,panels,400.,500.,650.,800.,outer]
    idx=[None]+list(range(len(parts)))+[None]*5
    assert len(idx)==len(edges)-1
    t,Dt=grid(n);blocks=len(idx);size=blocks*(n+1)
    A=lil_matrix((size,size),dtype=complex);b=np.zeros(size,complex);meta=[]
    for j,(lo,hi) in enumerate(zip(edges[:-1],edges[1:])):
        zl,zh=np.log([lo-rp,hi-rp]);z=(zl+zh)/2+(zh-zl)/2*t
        x=np.exp(z);rad=rp+x;Dz=Dt*2/(zh-zl);Dzz=Dz@Dz
        delta,dp,V=radial_coefficients(rad,a,mu,omega,m,lam)
        # Multiply the transformed ODE by x^2/Delta to keep its leading
        # coefficient one. Derivatives and interface matching remain in r.
        B=Dzz+((dp*x/delta-1)[:,None])*Dz+np.diag(V*x*x/delta)
        rhs=source_eval(rad,idx[j],parts,rp)*x*x/delta
        sl=slice(j*(n+1),(j+1)*(n+1));A[sl,sl]=B;b[sl]=rhs
        meta.append((rad,Dz/x[:,None],zl,zh))
    # At each shared boundary replace the two duplicate PDE rows with
    # continuity and radial-derivative continuity. J may jump there.
    for j in range(blocks-1):
        end=(j+1)*(n+1)-1;start=end+1
        A[end,:]=0;A[end,end]=1;A[end,start]=-1;b[end]=0
        A[start,:]=0
        A[start,j*(n+1):(j+1)*(n+1)]=meta[j][1][-1]
        A[start,(j+1)*(n+1):(j+2)*(n+1)]=-meta[j+1][1][0];b[start]=0
    hin,hd=horizon_series(edges[0],a,mu,omega,m,lam,order=4)
    up,ud,_=formal_boundary(a,mu,omega,m,lam,outer,order=4,dps=64,propagating=False)
    A[0,:]=0;A[0,0:n+1]=meta[0][1][0];A[0,0]-=hd/hin;b[0]=0
    A[-1,:]=0;A[-1,-n-1:]=meta[-1][1][-1];A[-1,-1]-=ud/up;b[-1]=0
    A=A.tocsr();rowscale=np.asarray(abs(A).max(axis=1).toarray()).ravel()
    scaled=A.multiply((1/rowscale)[:,None]).tocsr();bb=b/rowscale
    sol=spsolve(scaled,bb)
    zh=sol[0]/hin
    flux=mode_flux(omega,m,d['cloud_provenance']['temporal_omega'],d['cloud'],mu,a,0j,zh)['horizon']['orbital_energy']/mu**6
    # Independently evaluate differential residuals at non-collocation points.
    residuals=[];interfaces=[]
    for j,(rad,Dr,zl,zr) in enumerate(meta):
        values=sol[j*(n+1):(j+1)*(n+1)];coef=C.chebfit(t,values,n)
        test=np.cos(np.pi*(np.arange(n)+.5)/n);z=(zl+zr)/2+(zr-zl)/2*test
        x=np.exp(z);rr=rp+x;v=C.chebval(test,coef)
        vz=C.chebval(test,C.chebder(coef))*2/(zr-zl)
        vzz=C.chebval(test,C.chebder(coef,2))*(2/(zr-zl))**2
        delta,dp,V=radial_coefficients(rr,a,mu,omega,m,lam)
        J=source_eval(rr,idx[j],parts,rp)
        terms=np.array([delta/x**2*(vzz-vz),dp/x*vz,V*v,-J])
        residuals.append(float(np.linalg.norm(terms.sum(axis=0))/max(sum(np.linalg.norm(q) for q in terms),1e-300)))
    return dict(degree=n,unknowns=size,horizon_amplitude=[zh.real,zh.imag],horizon_flux_Li_units=float(flux),max_panel_offgrid_scaled_residual=max(residuals),linear_scaled_residual=float(np.linalg.norm(scaled@sol-bb)/np.linalg.norm(bb)))

def main():
    report={'method':'multidomain Chebyshev-Lobatto BVP in log(r-rplus)','scope':'frozen-source radial solver control; not a metric reconstruction replacement','modes':[],'limitations':['Same stored finite-resolution source and angular eigenvalue as production.','Same finite-end order-four boundary prescriptions; does not independently validate asymptotics.','Source interpolation fixed per original Gauss panel; spectral-order convergence is not source convergence.','Finite domain, no infinity compactification or fitted amplitude.']}
    for c,l,m in [(2,1,1),(2,3,1),(2,0,0)]:
        path=ROOT/f'outputs/li_fig9_10_20260928/r20_c{c}_l{l}_m{m}.json'
        raw=path.read_bytes();d=json.loads(raw);rp=1+np.sqrt(1-.88**2)
        panels,r,J,parts=source_parts(d,rp)
        green=LiOrder4Green(.88,.3,d['omega'],l,m)
        response=SampledResponse(green,panels,r,J,log_first=True,rtol=2e-12,atol=1e-14)
        refzh=response.horizon_coefficient
        refflux=mode_flux(d['omega'],m,d['cloud_provenance']['temporal_omega'],c,.3,.88,0j,refzh)['horizon']['orbital_energy']/.3**6
        entry={'cloud':c,'ell':l,'m':m,'input':str(path.relative_to(ROOT)),'input_sha256':hashlib.sha256(raw).hexdigest(),'stored_Gauss_flux':d['flux']['horizon'],'continuous_Green_flux':refflux,'continuous_Green_vs_stored_flux_relative':refflux/d['flux']['horizon']-1,'trials':[]}
        for n in [16,24,32,48,64,80,96]:
            start=time.monotonic();row=spectral(d,n,parts,panels)
            zh=complex(*row['horizon_amplitude']);row.update(complex_amplitude_relative_error_vs_Green=float(abs(zh/refzh-1)),flux_relative_error_vs_Green=float(row['horizon_flux_Li_units']/refflux-1),seconds=time.monotonic()-start)
            entry['trials'].append(row);print(c,l,m,n,row['flux_relative_error_vs_Green'],row['max_panel_offgrid_scaled_residual'],flush=True)
        report['modes'].append(entry)
        (OUT/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    report['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (OUT/'results.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
