"""Independent bound scalar00 Up/W kernel via high-precision In series.

At real frequency In and conjugate(In) span the radial equation. Match their
combination to an independently generated decaying asymptotic series at a
finite distant radius. This avoids the production Up ODE solver entirely.
"""
import json
from pathlib import Path
import mpmath as mp
import numpy as np
from paper_leaver_radial import LeaverIngoing
from environment_radial import RadialGreen


def decaying_logder(r,a,mu,w,separation,order):
    q=-mp.sqrt(mu*mu-w*w);p=(mu*mu-2*w*w)/q-1
    d2=[a**4,-4*a*a,4+2*a*a,-4,1];dd=[-2*a*a,4+2*a*a,-6,2]
    b=a*a*w;v=[b*b-a*a*separation,2*separation,2*w*b-separation-a*a*mu*mu,2*mu*mu,w*w-mu*mu]
    def entry(n,target):
        value=mp.mpf(0)
        for degree,c in enumerate(d2):
            for shift,factor in [(0,q*q),(1,2*q*(p-n)),(2,(p-n)*(p-n-1))]:
                if degree-n-shift==target:value+=c*factor
        for degree,c in enumerate(dd):
            if degree-n==target:value+=c*q
            if degree-n-1==target:value+=c*(p-n)
        for degree,c in enumerate(v):
            if degree-n==target:value+=c
        return value
    coefficients=[mp.mpf(1)]
    for n in range(1,order+1):
        coefficients.append(-sum(entry(j,3-n)*c for j,c in enumerate(coefficients))/entry(n,3-n))
    f=sum(c/r**n for n,c in enumerate(coefficients))
    fp=sum(-n*c/r**(n+1) for n,c in enumerate(coefficients))
    return q+p/r+fp/f,abs(coefficients[-1]/r**order/f)


def main():
    root=Path(__file__).resolve().parents[1];folder=root/'docs/environment_reproduction'
    old=json.loads((folder/'forced_mode_nr8_nt18_L18_mg-1_sl0_inner0.0005_outer320_log_h32.json').read_text())
    p=old['parameters'];a=p['metric']['a'];mu=p['alpha'];w=p['omega']
    green=RadialGreen(a,mu,w,0,0,rmax=1000.,offset=1e-4,rtol=1e-11)
    series=LeaverIngoing(a,mu,w,0,0,terms=60000,dps=90);rows=[]
    with mp.workdps(90):
        aa,mm,ww=[mp.mpf(str(x)) for x in (a,mu,w)];separation=series.angular+aa*aa*ww*ww
        for match,order in [(150,20),(200,24),(200,32),(250,40)]:
            ri,rip=series.state_mp(match);u0=mp.matrix([[ri,mp.conj(ri)],[rip,mp.conj(rip)]])
            logarithmic,last=decaying_logder(mp.mpf(match),aa,mm,ww,separation,order)
            amplitudes=mp.lu_solve(u0,mp.matrix([1,logarithmic]));samples=[]
            for r in [green.rp+.001,3.,10.,20.,40.]:
                inside,dinside=series.state_mp(r);up=amplitudes[0]*inside+amplitudes[1]*mp.conj(inside)
                dup=amplitudes[0]*dinside+amplitudes[1]*mp.conj(dinside)
                delta=(mp.mpf(str(r))-series.rp)*(mp.mpf(str(r))-series.rm)
                W=delta*(inside*dup-up*dinside)
                kernel=complex(up/W);reference=green.upsol.sol(r)[0]/green.w0
                samples.append(dict(r=r,relative_kernel_error=float(abs(kernel/reference-1))))
            row=dict(match=match,asymptotic_order=order,last_term_relative=float(last),samples=samples)
            rows.append(row);print(match,order,max(s['relative_kernel_error'] for s in samples),flush=True)
    result=dict(status='independent_bound_scalar00_Up_over_W_check',a=a,mu=mu,omega=w,rows=rows,
       limitations=['One dominant bound scalar00 channel; not all propagating modes',
         'Finite distant matching asymptotics need cutoff/order variation',
         'No production Up integrator, radial boundary helper or source used',
         'No flux fit or adjustable amplitude normalization'])
    (folder/'leaver_green_kernel_audit.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
