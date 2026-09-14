"""Actual 2023 ten-component spin-weighted tetrad matching diagnostic.

All equations are held out: no gauge coefficients are fitted here.
Includes the rho h_lplus_mplus component required for m=1.
"""
import argparse,json
from pathlib import Path
import numpy as np
from pybhpt.swsh import Yslm
from lorenz_ghp import KerrGHP
from lorenz_jet import Jet
from environment_source import kerr_metric

PAIRS=[(a,b) for a in range(4) for b in range(a,4)]
SPINS=(0,0,2,-2,1,-1,1,-1,0,0)
LABELS=('lp_lp','lm_lm','mp_mp','mm_mm','rho_lp_mp','rhobar_lp_mm','rhobar_lm_mp','rho_lm_mm','SigmaDelta_lp_lm','trace')


def transform(r,t,a):
    g=KerrGHP(r,t,a,order=2);z=Jet(0.,2);one=Jet(1.,2)
    lp=[(g.r*g.r+a*a)/g.delta,one,z,a/g.delta]
    lm=[-lp[0],one,z,-lp[3]]
    mp=[1j*a*g.theta.sin(),z,one,1j/g.theta.sin()]
    mm=[-mp[0],z,one,-mp[3]]
    rho=g.zeta.conjugate();rhobar=g.zeta
    rows=[]
    for u,v,factor in ((lp,lp,1),(lm,lm,1),(mp,mp,1),(mm,mm,1),
        (lp,mp,rho),(lp,mm,rhobar),(lm,mp,rhobar),(lm,mm,rho),(lp,lm,g.sigma*g.delta)):
        rows.append([factor*(u[i]*v[j]+(u[j]*v[i] if i!=j else 0)) for i,j in PAIRS])
    rows.append([g.inv[i][j]*(2 if i!=j else 1) for i,j in PAIRS])
    return np.array([[[v.value for v in row] for row in rows],
                     [[v.derivative(0).value for v in row] for row in rows]])


def main():
    p=argparse.ArgumentParser();p.add_argument('--ellmax',type=int,default=8)
    p.add_argument('--testmax',type=int,default=6);p.add_argument('--quadrature',type=int,default=24)
    p.add_argument('--resume',action='store_true')
    args=p.parse_args()
    from environment_angular_diagnostic import install_dense_angular_diagnostic
    from environment_trace_variation import install_analytic_kappa_diagnostic
    install_dense_angular_diagnostic();install_analytic_kappa_diagnostic()
    Jet.coefficient_dtype=np.clongdouble
    from lorenz_metric import nonstatic_metric
    r0=20.;a=.8771530275949366;m=1;eps=5e-5
    x,w=np.polynomial.legendre.leggauss(args.quadrature);theta=np.arccos(x)
    degrees=list(range(1,args.testmax+1))
    angular=np.array([[Yslm(s,j,m,theta) if j>=max(abs(s),abs(m)) else np.zeros_like(theta) for s in SPINS] for j in degrees])
    weights=2*np.pi*angular*w
    T=np.array([transform(r0,t,a) for t in theta])
    metric=kerr_metric(r0,np.pi/2,a);op=1/(r0**1.5+a)
    ut=1/np.sqrt(-metric[0,0]-2*op*metric[0,3]-op*op*metric[3,3]);ucov=metric@np.array([ut,0,0,op*ut])
    coordinate_target=-8*(np.outer(ucov,ucov)+metric/2)/(ut*(r0*r0-2*r0+a*a))
    projected_target=transform(r0,np.pi/2,a)[0]@np.array([coordinate_target[i,j] for i,j in PAIRS])
    target=np.array([[2*np.pi*Yslm(s,j,m,np.pi/2)*projected_target[c] if j>=max(abs(s),abs(m)) else 0 for c,s in enumerate(SPINS)] for j in degrees])
    sums=np.zeros((2,len(degrees),10),complex);rows=[]
    folder=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    out=folder/f'nonstatic_paper_projection_r20_m1_q{args.quadrature}_j{args.testmax}.json'
    scale=np.maximum(np.max(abs(target),axis=0),1e-12)
    first=1
    if args.resume:
        old=json.loads(out.read_text())
        required=dict(r0=r0,a=a,m=m,quadrature=args.quadrature,epsilon=eps,test_degrees=degrees,
                      tetrad_convention='2306.16459 Eq tetrad2; lm past-directed')
        if any(old.get(k)!=v for k,v in required.items()):raise ValueError('Resume configuration mismatch')
        def dec(v):
            v=np.asarray(v);return v[...,0]+1j*v[...,1]
        np.testing.assert_allclose(dec(old['target']),target,rtol=2e-13,atol=2e-12)
        rows=old['rows'];last=rows[-1]
        sums[0]=dec(last['value_jump']);sums[1]=dec(last['derivative_error'])+target
        first=last['metric_ellmax']+1
    for L in range(first,args.ellmax+1):
        delta=np.zeros_like(sums)
        for sign in (-1,1):
            values=[]
            for t in theta:
                _,h=nonstatic_metric(r0+sign*eps,t,r0,a,L,m,order=8)
                jet=np.array([[getattr(h[i][j],'value') for i,j in PAIRS],
                    [h[i][j].derivative(0).value for i,j in PAIRS],
                    [h[i][j].derivative(0).derivative(0).value for i,j in PAIRS]])
                dr=-sign*eps;values.append([jet[0]+dr*jet[1]+.5*dr*dr*jet[2],jet[1]+dr*jet[2]])
            values=np.array(values,dtype=complex)
            h0=np.einsum('kfc,kc->kf',T[:,0],values[:,0])
            h1=np.einsum('kfc,kc->kf',T[:,0],values[:,1])+np.einsum('kfc,kc->kf',T[:,1],values[:,0])
            delta[0]+=sign*np.einsum('jfk,kf->jf',weights,h0)
            delta[1]+=sign*np.einsum('jfk,kf->jf',weights,h1)
        sums+=delta
        rows.append(dict(metric_ellmax=L,value_jump=enc(sums[0]),derivative_error=enc(sums[1]-target),
            value_jump_scaled=float(np.max(abs(sums[0])/scale)),derivative_error_scaled=float(np.max(abs(sums[1]-target)/scale)),
            m1_extra_component_j1=dict(value_jump=enc(sums[0,0,4]),derivative_error=enc(sums[1,0,4]-target[0,4]))))
        result=dict(status='held_out_paper_tetrad_matching_diagnostic',tetrad_convention='2306.16459 Eq tetrad2; lm past-directed',r0=r0,a=a,m=m,quadrature=args.quadrature,
            epsilon=eps,test_degrees=degrees,components=LABELS,spins=SPINS,target=enc(target),component_scales=np.asarray(scale,float).tolist(),rows=rows,
            limitations=['Finite separated-mode sum; tests are not independent of truncation',
                         'Scaled residual uses largest target in each component, not componentwise relative error',
                         'No coefficients fitted; scalar flux unchanged'])
        tmp=out.with_suffix('.tmp');tmp.write_text(json.dumps(result,indent=2)+'\n');tmp.replace(out)
        print(L,rows[-1]['value_jump_scaled'],rows[-1]['derivative_error_scaled'],rows[-1]['m1_extra_component_j1'],flush=True)


def enc(a):
    a=np.asarray(a,dtype=complex);return np.stack([a.real,a.imag],axis=-1).tolist()


if __name__=='__main__':main()
