"""Independent compact-support quadrature of N^dag=-N (bilinear, no conjugate)."""
import json
from pathlib import Path
import numpy as np
from lorenz_ghp import KerrGHP
from lorenz_corrector import corrector


def check(n):
    x,w=np.polynomial.legendre.leggauss(n)
    rr,tt=6+2*x,1+.5*x
    total,scale=0j,0.
    A=np.array([[1,.2,.1,.3],[.2,.5,.4,.1],[.1,.4,.2,.6],[.3,.1,.6,.7]])
    B=np.array([[.4,.3,.2,.1],[.3,.8,.1,.4],[.2,.1,.5,.2],[.1,.4,.2,.6]])
    for i,r in enumerate(rr):
        for j,t in enumerate(tt):
            g=KerrGHP(r,t,.6,omega=.13,m=2,order=4)
            f=(g.r-4)**3*(8-g.r)**3*(g.theta-.5)**3*(1.5-g.theta)**3
            u=[[f*A[a,b] for b in range(4)] for a in range(4)]
            v=[[f*(1+.1*g.r+.2j*g.theta)*B[a,b] for b in range(4)] for a in range(4)]
            nu=np.array([[z.value for z in row] for row in corrector(g,u)])
            g.omega,g.m=-.13,-2
            nv=np.array([[z.value for z in row] for row in corrector(g,v)])
            inv=np.array([[z.value for z in row] for row in g.inv])
            uv=np.array([[z.value for z in row] for row in u])
            vv=np.array([[z.value for z in row] for row in v])
            left=np.sum((inv@vv@inv)*nu)
            right=np.sum((inv@uv@inv)*nv)
            weight=w[i]*w[j]*(g.sigma*g.theta.sin()).value
            total+=weight*(left+right)
            scale+=abs(weight)*(abs(left)+abs(right))
    return dict(n=n,relative_defect=abs(total)/scale,integral=[total.real,total.imag],scale=scale)


if __name__=='__main__':
    reports=[]
    for n in (8,12):
        result=check(n)
        reports.append(result)
        print(result,flush=True)
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/corrector_adjoint.json'
    out.write_text(json.dumps(reports,indent=2)+'\n')
