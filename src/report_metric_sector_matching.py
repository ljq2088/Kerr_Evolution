"""Resolve high-L projected jump increments into spin reconstruction sectors."""
import json
import hashlib
from pathlib import Path
import numpy as np
from scipy.special import lpmv
from lorenz_metric import spin2_metric,spin1_metric,spin0_metric
from environment_trace_variation import install_analytic_kappa_diagnostic


def main():
    install_analytic_kappa_diagnostic()
    directory=Path(__file__).resolve().parents[1]/'docs/environment_reproduction'
    reference=directory/'metric_tensor_matching_q18_smooth_m1_r10_a0.877153027595_analytickappa.json'
    ref=json.loads(reference.read_text())
    assert [c['ellmax'] for c in ref['cases']]==list(range(1,19))
    a=ref['a'];r0=ref['r0'];m=ref['m'];epsilon=ref['epsilon']
    x,w=np.polynomial.legendre.leggauss(ref['quadrature'])
    indices=[tuple(v) for v in ref['components']]
    power=np.array([int(i==2)+int(j==2) for i,j in indices])
    weights=np.array([w*lpmv(abs(m),ell,x) for ell in (abs(m),abs(m)+1)])[:,:,None]*(1-x*x)[None,:,None]**(power[None,None,:]/2)
    def encode(v):return np.stack((v.real,v.imag),axis=-1).tolist()
    def decode(v):
        v=np.asarray(v);return v[...,0]+1j*v[...,1]
    cases=[]
    for ell in (17,18):
        sectors=[]
        for name,function in (('spin2',spin2_metric),('spin1',spin1_metric),('spin0',spin0_metric)):
            sides=[]
            for sign in (-1,1):
                r=r0+sign*epsilon;acc=np.zeros((2,3,10),complex)
                for node,theta in enumerate(np.arccos(x)):
                    kwargs={'full_current':True} if name=='spin1' else {}
                    _,h=function(r,theta,r0,a=a,ell=ell,m=m,order=8,**kwargs)
                    values=np.array([[h[i][j].value for i,j in indices],
                        [h[i][j].derivative(0).value for i,j in indices],
                        [h[i][j].derivative(0).derivative(0).value for i,j in indices]])
                    acc+=weights[:,node,None,:]*values[None,:,:]
                d=-sign*epsilon
                sides.append((acc[:,0]+d*acc[:,1]+d*d*acc[:,2]/2,acc[:,1]+d*acc[:,2]))
            v=sides[1][0]-sides[0][0];dv=sides[1][1]-sides[0][1]
            sectors.append(dict(name=name,value_jump=encode(v),derivative_jump=encode(dv),maximum_value_jump=float(np.max(abs(v))),maximum_derivative_jump=float(np.max(abs(dv)))))
            print(ell,name,sectors[-1]['maximum_value_jump'],sectors[-1]['maximum_derivative_jump'],flush=True)
        summed=[sum(decode(s[key]) for s in sectors) for key in ('value_jump','derivative_jump')]
        previous,current=[ref['cases'][k]['limits'][-1] for k in (ell-2,ell-1)]
        changes=[decode(current[key])-decode(previous[key]) for key in ('value_jump','derivative_jump')]
        cases.append(dict(ell=ell,sectors=sectors,total_value_jump=encode(summed[0]),total_derivative_jump=encode(summed[1]),
            original_sequence_value_difference=float(np.max(abs(summed[0]-changes[0]))),
            original_sequence_derivative_difference=float(np.max(abs(summed[1]-changes[1])))))
        out=directory/'rp10_metric_high_ell_sector_audit.json';tmp=out.with_suffix('.tmp')
        tmp.write_text(json.dumps(dict(status='complete' if ell==18 else 'in_progress',
            reference=reference.name,reference_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),
            a=a,r0=r0,m=m,quadrature=ref['quadrature'],epsilon=epsilon,kappa_backend='analytic',cases=cases,
            limitation='Individual sectors need not satisfy physical matching. These are L17 and L18 increments, not the complete metric or a scalar-flux error estimate.'),indent=2)+'\n');tmp.replace(out)


if __name__=='__main__':main()
