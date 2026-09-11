"""Local horizon integral from the two Frobenius phases, real nonzero gamma.

Fits the radial Green integrand, not a desired flux. Order/window changes and
held-out samples are diagnostics; no rigorous extrapolation bound is claimed.
"""
import numpy as np


def horizon_tail(r,integrand,rp,gamma,cutoff,order=1,window=10.):
    if not np.isreal(gamma) or abs(gamma)<1e-5 or cutoff<=0 or order<0 or window<=1:
        raise ValueError('Nonzero real horizon frequency and a positive fitting interval required')
    r=np.asarray(r);values=np.asarray(integrand,complex)
    if r.ndim!=1 or values.shape!=r.shape or np.any(np.diff(r)<=0):
        raise ValueError('Sorted radii and corresponding integrand required')
    selected=(r-rp>=cutoff*(1-1e-10))&(r-rp<=window*cutoff)
    x=(r[selected]-rp)/cutoff;y=values[selected]
    powers=np.array([n+phase for phase in (0.,-2j*gamma) for n in range(order+1)])
    if len(x)<len(powers)+2:
        raise ValueError('Insufficient near-horizon samples for this order and window')
    design=np.exp(np.log(x)[:,None]*powers[None,:])
    norms=np.linalg.norm(design,axis=0)
    scaled=design/norms
    coefficient,_,rank,singular=np.linalg.lstsq(scaled,y,rcond=1e-12)
    if rank<len(powers):raise ValueError('Horizon fit is rank deficient')
    coefficient/=norms
    correction=cutoff*np.sum(coefficient/(powers+1))
    residual=np.max(abs(design@coefficient-y))
    holdout=np.arange(len(x))%3==1
    holdout_residual=None
    if np.sum(~holdout)>=len(powers)+1:
        training,_,training_rank,_=np.linalg.lstsq(scaled[~holdout],y[~holdout],rcond=1e-12)
        if training_rank==len(powers):
            holdout_residual=float(np.max(abs(scaled[holdout]@training-y[holdout])))
    return correction,dict(order=order,window=window,samples=int(len(x)),rank=int(rank),
        column_normalized_condition=float(singular[0]/singular[-1]),
        maximum_integrand_fit_residual=float(residual),held_out_integrand_residual=holdout_residual,
        maximum_sample_integrand=float(np.max(abs(y))),
        method='two horizon Frobenius phases; analytic integration from zero to cutoff')
