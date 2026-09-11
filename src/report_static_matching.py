"""Local static junction diagnostic; no horizon/infinity boundary claim.

Fit homogeneous scalar Hessian and completion jumps to projected particle
conditions. Coefficients describe outside-minus-inside local solutions only.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import eval_legendre
from environment_source import kerr_metric
from lorenz_ghp import KerrGHP
from lorenz_jet import Jet
from lorenz_tensor import vector_covariant_derivative
from lorenz_static_gauge import _angular_jet, static_trace_metric
from lorenz_static_spin2 import sourced_static_spin2
from lorenz_completion import completion_metric


def scalar_hessian(r,theta,a,ell,datum,order=6):
    """Hessian of Box kappa=0 with local (R,R')=(1,0) or (0,1)."""
    g=KerrGHP(r,theta,a,order=order)
    radial=Jet(float(datum==0),order)
    radial.c[1,0]=float(datum==1)
    for n in range(order-1):
        rhs=(ell*(ell+1)*radial-2*(g.r-1)*radial.derivative(0))/g.delta
        radial.c[n+2,0]=rhs.c[n,0]/((n+1)*(n+2))
    field=radial*_angular_jet(g,ell)
    return g,vector_covariant_derivative(g,[g.partial(field,j) for j in range(4)])


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--a',type=float,default=.6)
    p.add_argument('--r0',type=float,default=6.)
    p.add_argument('--ellmax',type=int,default=4)
    p.add_argument('--quadrature',type=int,default=10)
    p.add_argument('--testmax',type=int,default=2)
    p.add_argument('--epsilon',type=float,default=5e-5)
    p.add_argument('--continuity-only',action='store_true',help='Fit values and conserved charges; hold out all derivative jumps')
    p.add_argument('--free-ellmax',type=int,default=None,help='Largest free scalar degree; default ellmax+2')
    p.add_argument('--reuse-samples',action='store_true',help='Reuse matching vacuum jets at the same angular nodes and separation')
    args=p.parse_args()
    if args.ellmax<2 or args.testmax<0 or args.quadrature<args.testmax+3:
        raise ValueError('Require ellmax>=2 and adequate angular quadrature')
    x,w=np.polynomial.legendre.leggauss(args.quadrature)
    pairs=[(i,j) for i in range(4) for j in range(i,4)]
    # Include two held-out angular degrees, excluded from the least-squares fit.
    degrees=np.arange(args.testmax+3)
    weights=np.array([w*eval_legendre(j,x) for j in degrees])[:,:,None]
    weights=weights*(np.sqrt(1-x*x)[:,None]**np.array([int(i==2)+int(j==2) for i,j in pairs]))[None,:,:]
    shape=(2,len(degrees),10)
    particular=np.zeros(shape,complex)
    def values(h,shift=0.):
        v=np.array([h[i][j].value for i,j in pairs])
        d=np.array([h[i][j].derivative(0).value for i,j in pairs])
        dd=np.array([h[i][j].derivative(0).derivative(0).value for i,j in pairs])
        return np.array([v+shift*d+shift*shift*dd/2,d+shift*dd])
    root=Path(__file__).resolve().parents[1]
    sample_path=root/'outputs'/f'static_samples_a{args.a:g}_r{args.r0:g}_q{args.quadrature}_eps{args.epsilon:g}.npz'
    sample_metadata=dict(version=1,a=args.a,r0=args.r0,quadrature=args.quadrature,epsilon=args.epsilon,
                         order=8,circular_isometry_average=True)
    sampled_particular={}
    sampled_basis={}
    if args.reuse_samples and sample_path.exists():
        with np.load(sample_path,allow_pickle=False) as saved:
            if json.loads(str(saved['metadata']))!=sample_metadata:
                raise ValueError('Saved samples have different physics or discretization')
            sampled_particular={int(k[2:]):saved[k] for k in saved.files if k.startswith('p_')}
            sampled_basis={k[2:]:saved[k] for k in saved.files if k.startswith('b_')}
    def save_samples():
        temporary=sample_path.with_suffix('.tmp')
        with temporary.open('wb') as stream:
            np.savez_compressed(stream,metadata=json.dumps(sample_metadata),
                **{f'p_{ell}':v for ell,v in sampled_particular.items()},
                **{f'b_{label}':v for label,v in sampled_basis.items()})
        temporary.replace(sample_path)
    for ell in range(2,args.ellmax+1):
        if ell in sampled_particular:
            particular+=np.einsum('dkc,jkc->djc',sampled_particular[ell],weights)
            print(f'Reused particular ell={ell}',flush=True)
            continue
        nodal=np.zeros((2,len(x),10),complex)
        for sign in (-1,1):
            r=args.r0+sign*args.epsilon
            for k,t in enumerate(np.arccos(x)):
                _,h=sourced_static_spin2(r,t,r0=args.r0,a=args.a,ell=ell,order=8)
                v=values(h,-sign*args.epsilon)
                if ell%2==0:
                    _,h,_=static_trace_metric(r,t,r0=args.r0,a=args.a,ell=ell,order=8)
                    v+=values(h,-sign*args.epsilon)
                nodal[:,k]+=sign*v
        sampled_particular[ell]=nodal
        save_samples()
        particular+=np.einsum('dkc,jkc->djc',nodal,weights)
        print(f'Completed particular ell={ell}',flush=True)
    free_ellmax=args.ellmax+2 if args.free_ellmax is None else args.free_ellmax
    labels=list('BCDEFG')+[f'kappa_{ell}_{d}' for ell in range(2,free_ellmax+1,2) for d in (0,1)]
    basis=np.zeros((len(labels),)+shape,complex)
    for n,label in enumerate(labels):
        if label in sampled_basis:
            basis[n]=np.einsum('dkc,jkc->djc',sampled_basis[label],weights)
            continue
        nodal=np.zeros((2,len(x),10),complex)
        for k,t in enumerate(np.arccos(x)):
            if len(label)==1:
                _,h,_=completion_metric(args.r0,t,a=args.a,mode=label,order=6,reference_radius=args.r0)
            else:
                _,ell,d=label.split('_')
                _,h=scalar_hessian(args.r0,t,args.a,int(ell),int(d))
            nodal[:,k]=values(h)
        sampled_basis[label]=nodal
        save_samples()
        basis[n]=np.einsum('dkc,jkc->djc',nodal,weights)
    metric=kerr_metric(args.r0,np.pi/2,args.a)
    op=1/(args.r0**1.5+args.a)
    ut=1/np.sqrt(-metric[0,0]-2*op*metric[0,3]-op*op*metric[3,3])
    u=metric@np.array([ut,0,0,op*ut])
    jump=-8*(np.outer(u,u)+metric/2)/(ut*(args.r0**2-2*args.r0+args.a**2))
    target=np.zeros(shape)
    target[1]=np.array([[eval_legendre(j,0)*jump[i,k] for i,k in pairs] for j in degrees])
    # Scale coordinate components and radial derivatives dimensionlessly.
    lengths=np.array([1,1,args.r0,args.r0])
    scale=np.array([1/(lengths[i]*lengths[j]) for i,j in pairs])[None,None,:]*np.array([1,args.r0])[:,None,None]
    used=(slice(0,1) if args.continuity_only else slice(None),slice(0,args.testmax+1),slice(None))
    matrix=np.moveaxis(basis.real*scale,0,-1)[used].reshape(-1,len(labels))
    rhs=((target-particular.real)*scale)[used].ravel()
    if args.continuity_only:
        # Eliminate the known charges exactly; no least-squares relaxation.
        prescribed={labels.index('E'):float(-u[0]),labels.index('G'):float(u[3]+args.a*u[0])}
    else:
        prescribed={}
    for index,value in prescribed.items():
        rhs-=matrix[:,index]*value
    norms=np.linalg.norm(matrix,axis=0)
    active=norms>1e-9*np.max(norms)
    for index in prescribed:
        active[index]=False
    coefficient=np.zeros(len(labels))
    for index,value in prescribed.items():
        coefficient[index]=value
    fitted,_,rank,singular=np.linalg.lstsq(matrix[:,active]/norms[active],rhs,rcond=1e-10)
    coefficient[active]=fitted/norms[active]
    residual=particular.real+np.einsum('n,ndjc->djc',coefficient,basis.real)-target
    result=dict(status='local_junction_diagnostic_not_boundary_matched',parameters=vars(args),
        circular_isometry_average=True,
        derivative_conditions_used_in_fit=not args.continuity_only,
        charge_conditions='exact_elimination' if args.continuity_only else 'not_prescribed',
        residual_scale='value/(L_i L_j), r0*radial_derivative/(L_i L_j); L=(1,1,r0,r0)',
        components=pairs,test_degrees=degrees.tolist(),fit_testmax=args.testmax,
        basis_labels=labels,coefficients=coefficient.tolist(),rank=int(rank),
        unconstrained_columns=[label for n,(label,keep) in enumerate(zip(labels,active)) if not keep and n not in prescribed],
        singular_values_column_normalized=singular.tolist(),
        expected_charge_jumps=dict(E=float(-u[0]),G=float(u[3]+args.a*u[0])),
        fitted_completion_jumps={label:float(coefficient[n]) for n,label in enumerate(labels[:6])},
        max_imaginary_input=float(max(np.max(abs(particular.imag)),np.max(abs(basis.imag)))),
        max_scaled_fit_residual=float(np.max(abs((residual*scale)[used]))),
        max_scaled_derivative_residual_low_degrees=float(np.max(abs((residual*scale)[1,:args.testmax+1]))),
        max_scaled_holdout_residual=float(np.max(abs((residual*scale)[:,args.testmax+1:]))),
        scaled_value_residual_by_degree=np.max(abs((residual*scale)[0]),axis=-1).tolist(),
        scaled_derivative_residual_by_degree=np.max(abs((residual*scale)[1]),axis=-1).tolist(),
        particular_jump=particular.real.tolist(),target=target.tolist(),residual=residual.tolist())
    suffix='_continuity' if args.continuity_only else ''
    if args.free_ellmax is not None:
        suffix+=f'_free{args.free_ellmax}'
    if args.continuity_only:
        suffix+='_exactcharges'
    if args.epsilon!=5e-5:
        suffix+=f'_eps{args.epsilon:g}'
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'/f'static_matching_a{args.a:g}_r{args.r0:g}_L{args.ellmax}_q{args.quadrature}_j{args.testmax}{suffix}.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('rank','max_scaled_fit_residual','max_scaled_derivative_residual_low_degrees','max_scaled_holdout_residual','expected_charge_jumps','fitted_completion_jumps')},indent=2),flush=True)


if __name__=='__main__':
    main()
