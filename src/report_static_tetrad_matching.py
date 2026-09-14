"""Independent static junction check in the paper's spin-weighted tetrad basis.

Uses already sampled coordinate vacuum jets; fits continuity and exact charges,
and holds all derivative conditions out. No physical boundary claim.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.special import lpmv
from lorenz_ghp import KerrGHP
from lorenz_jet import Jet
from environment_source import kerr_metric


def tetrad_weights(r,theta,a):
    g=KerrGHP(r,theta,a,order=2)
    z=Jet(0.,2)
    lp=[(g.r*g.r+a*a)/g.delta,Jet(1.,2),z,a/g.delta]
    lm=[-lp[0],lp[1],z,-lp[3]]
    mp=[1j*a*g.theta.sin(),z,Jet(1.,2),1j/g.theta.sin()]
    pairs=[(i,j) for i in range(4) for j in range(i,4)]
    rows=[]
    for u,v,factor in ((lp,lp,1),(mp,mp,1),(lp,mp,g.zeta.conjugate()),(lp,lm,g.sigma*g.delta)):
        rows.append([factor*(u[i]*v[j]+(u[j]*v[i] if i!=j else 0)) for i,j in pairs])
    rows.append([g.inv[i][j]*(2 if i!=j else 1) for i,j in pairs])
    return np.array([[[v.value for v in row] for row in rows],
                     [[v.derivative(0).value for v in row] for row in rows]])


def spin_harmonic(ell,spin,x):
    if ell<spin:
        return np.zeros_like(x)
    norm=math.sqrt((2*ell+1)/(4*np.pi)*math.factorial(ell-spin)/math.factorial(ell+spin))
    return norm*lpmv(spin,ell,x)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--samples',required=True)
    p.add_argument('--ellmax',type=int,required=True)
    p.add_argument('--testmax',type=int,required=True)
    p.add_argument('--free-ellmax',type=int,required=True)
    p.add_argument('--paper-conditions',action='store_true',help='Fit only the independent continuity conditions listed in 2306.16459')
    args=p.parse_args()
    if args.testmax<2 or args.free_ellmax<2 or args.free_ellmax%2:
        raise ValueError('testmax>=2 and an even free-ellmax>=2 required')
    with np.load(args.samples,allow_pickle=False) as f:
        meta=json.loads(str(f['metadata']))
        part=sum(f[f'p_{ell}'] for ell in range(2,args.ellmax+1))
        labels=list('BCDEFG')+[f'kappa_{ell}_{datum}' for ell in range(2,args.free_ellmax+1,2) for datum in (0,1)]
        base=np.array([f[f'b_{label}'] for label in labels])
    r,a=meta['r0'],meta['a']
    x,w=np.polynomial.legendre.leggauss(meta['quadrature'])
    transform=np.array([tetrad_weights(r,theta,a) for theta in np.arccos(x)])
    def project(v):
        tv=np.array([np.einsum('kfc,kc->kf',transform[:,0],v[0]),
                     np.einsum('kfc,kc->kf',transform[:,0],v[1])+np.einsum('kfc,kc->kf',transform[:,1],v[0])])
        return np.einsum('dkf,jkf->djf',tv,angular)
    degrees=np.arange(args.testmax+3)
    spins=(0,2,1,0,0)
    angular=np.array([[spin_harmonic(j,s,x)*w for s in spins] for j in degrees]).transpose(0,2,1)
    particular=project(part)
    basis=np.array([project(v) for v in base])
    metric=kerr_metric(r,np.pi/2,a)
    omega=1/(r**1.5+a)
    ut=1/np.sqrt(-metric[0,0]-2*omega*metric[0,3]-omega**2*metric[3,3])
    ucov=metric@np.array([ut,0,0,omega*ut])
    jump=-8*(np.outer(ucov,ucov)+metric/2)/(ut*(r*r-2*r+a*a))
    pairs=[(i,j) for i in range(4) for j in range(i,4)]
    equatorial=tetrad_weights(r,np.pi/2,a)[0]@np.array([jump[i,j] for i,j in pairs])
    target=np.zeros_like(particular)
    target[1]=[[spin_harmonic(j,s,0.)*equatorial[k] for k,s in enumerate(spins)] for j in degrees]
    scale=np.array([1,1/r**2,1/r**2,1/(r*r*(r*r-2*r+a*a)),1])
    prescribed={labels.index('E'):float(-ucov[0]),labels.index('G'):float(ucov[3]+a*ucov[0])}
    known=sum(basis[k]*v for k,v in prescribed.items())
    free=[k for k in range(len(labels)) if k not in prescribed]
    selected=np.zeros(particular.shape[1:],bool)
    if args.paper_conditions:
        selected[0,0]=selected[0,3]=selected[1,2]=selected[2,2]=True
        for ell in range(2,args.free_ellmax+1,2):
            if ell>args.testmax:
                raise ValueError('Paper conditions need testmax >= free-ellmax')
            selected[ell,0]=selected[ell,1]=True
    else:
        selected[:args.testmax+1]=True
    rhs=(-(particular+known)[0]*scale)[selected]
    matrix=np.moveaxis(basis[free,0]*scale,0,-1)[selected]
    matrix=np.vstack((matrix.real,matrix.imag))
    rhs=np.concatenate((rhs.real,rhs.imag))
    norms=np.linalg.norm(matrix,axis=0)
    if np.any(norms<1e-14):
        raise ValueError('A chosen free mode is not constrained')
    fitted,_,rank,singular=np.linalg.lstsq(matrix/norms,rhs,rcond=1e-12)
    if rank!=len(free):
        raise ValueError(f'Underdetermined matching: rank {rank}, unknowns {len(free)}')
    coefficient=np.zeros(len(labels))
    coefficient[free]=fitted/norms
    for k,v in prescribed.items():coefficient[k]=v
    residual=particular+np.einsum('n,ndjf->djf',coefficient,basis)-target
    scaled=abs(residual*scale[None,None,:]*np.array([1,r])[:,None,None])
    encode=lambda v:np.stack((v.real,v.imag),axis=-1).tolist()
    result=dict(status='static_tetrad_local_diagnostic_not_boundary_matched',parameters=vars(args),sample_metadata=meta,
        components=['l+l+','m+m+','rho_l+m+','Sigma_Delta_l+l-','trace'],spins=spins,
        basis_labels=labels,coefficients=coefficient.tolist(),rank=int(rank),
        singular_values=singular.tolist(),charge_conditions='exact_elimination',derivative_conditions_used_in_fit=False,
        selected_continuity_conditions=selected.tolist(),
        maximum_selected_continuity_residual=float(np.max(scaled[0][selected])),
        scaled_value_residual_by_degree=np.max(scaled[0],axis=-1).tolist(),
        scaled_derivative_residual_by_degree=np.max(scaled[1],axis=-1).tolist(),
        particular=encode(particular),target=encode(target),residual=encode(residual),
        projected_basis=encode(basis))
    root=Path(__file__).resolve().parents[1]
    suffix='_paper' if args.paper_conditions else ''
    if meta['epsilon']!=5e-5:
        suffix+=f"_eps{meta['epsilon']:g}"
    out=root/'docs/environment_reproduction'/f'static_tetrad_a{a:g}_r{r:g}_L{args.ellmax}_q{len(x)}_j{args.testmax}_free{args.free_ellmax}{suffix}.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(out.name,flush=True)
    print('value',result['scaled_value_residual_by_degree'],flush=True)
    print('derivative',result['scaled_derivative_residual_by_degree'],flush=True)


if __name__=='__main__':
    main()
