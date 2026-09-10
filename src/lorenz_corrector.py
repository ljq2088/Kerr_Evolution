"""AAB corrector N, transcription of 2406.12510v3 Appendix B.

Acts on smooth test tensors. Distributional particle sources are evaluated
through its formal adjoint N^dag=-N, not by sampling a delta distribution.
"""
from lorenz_jet import Jet


class Weighted:
    def __init__(self,f,p=0,q=0):
        self.f,self.p,self.q=f,p,q
    def __add__(self,x):
        if not isinstance(x,Weighted):
            if x==0:
                return self
            x=Weighted(self.f.lift(x))
        if (self.p,self.q)!=(x.p,x.q):
            raise ValueError(f'GHP sum mismatch: {(self.p,self.q)} vs {(x.p,x.q)}')
        return Weighted(self.f+x.f,self.p,self.q)
    __radd__=__add__
    def __neg__(self):
        return Weighted(-self.f,self.p,self.q)
    def __sub__(self,x):
        return self+-x
    def __rsub__(self,x):
        return -self+x
    def __mul__(self,x):
        if isinstance(x,Weighted):
            return Weighted(self.f*x.f,self.p+x.p,self.q+x.q)
        return Weighted(self.f*x,self.p,self.q)
    __rmul__=__mul__
    def __truediv__(self,x):
        return self*(1/x)
    def conjugate(self):
        return Weighted(self.f.conjugate(),self.q,self.p)


def tensor_covariant_derivative(g,tensor):
    return [[[g.partial(tensor[a][b],c)
              -sum(g.gamma[k][c][a]*tensor[k][b]+g.gamma[k][c][b]*tensor[a][k] for k in range(4))
              for b in range(4)] for a in range(4)] for c in range(4)]


def corrector(g,tensor):
    """Return covariant BL N[tensor]; no 8pi factor is absorbed here."""
    weights=[(1,1),(-1,-1),(1,-1),(-1,1)]
    T={}
    for A in range(4):
        for B in range(4):
            T[A,B]=Weighted(sum(g.tetrad[A][a]*g.tetrad[B][b]*tensor[a][b]
                               for a in range(4) for b in range(4)),
                            weights[A][0]+weights[B][0],weights[A][1]+weights[B][1])
    nabla=tensor_covariant_derivative(g,tensor)
    div=[sum(g.inv[a][c]*nabla[c][a][b] for a in range(4) for c in range(4)) for b in range(4)]
    divt=[Weighted(sum(g.tetrad[A][a]*div[a] for a in range(4)),*weights[A]) for A in range(4)]
    z=Weighted(g.zeta**4)
    psi,psib=Weighted(g.psi2),Weighted(g.psi2.conjugate())
    def components(prime=False):
        l,n,m,b=(1,0,3,2) if prime else (0,1,2,3)
        names=('rhop','rho','taup','tau') if prime else ('rho','rhop','tau','taup')
        rw=[(-1,-1),(1,1),(-1,1),(1,-1)] if prime else [(1,1),(-1,-1),(1,-1),(-1,1)]
        rho,rhop,tau,taup=[Weighted(g.sc[name],*w) for name,w in zip(names,rw)]
        rb,rpb,tb,tpb=[x.conjugate() for x in (rho,rhop,tau,taup)]
        kinds=('thornp','thorn','ethp','eth') if prime else ('thorn','thornp','eth','ethp')
        def derivative(x,kind):
            val,p,q=g.derivative(x.f,x.p,x.q,kind)
            return Weighted(val,p,q)
        d,dp,e,ep=[lambda x,k=k:derivative(x,k) for k in kinds]
        trace=2*(T[m,b]-T[l,n])
        tilde=2*(T[m,b]+T[l,n])
        L,N,M,B=divt[l],divt[n],divt[m],divt[b]
        # Each bracket acts on everything to its right, including zeta^4.
        u=z*T[l,l]
        ll=(3*(ep(e(u)+(4*tau-tpb)*u)-taup*(e(u)+(4*tau-tpb)*u))
            +d(dp(u)+(4*rhop-rpb)*u)-rho*(dp(u)+(4*rhop-rpb)*u)
            -(2*rho*rhop+9*psi)*u)/9
        ll=ll-d(z*d(trace))/54+rho*(2*rho-rb)*z*tilde/9
        u=z*T[l,b]
        ll=ll+2*(tau*d(u)-rho*(e(u)+4*tau*u))/9
        u=z*T[l,m]
        inner=2*d(u)+(6*rho-3*rb)*u
        ll=ll-(2*ep(inner)-tb*inner-(8*rho*taup+5*rb*tb)*u)/9+2*d(z*L)/27

        u=z*T[m,m]
        mm=-(3*(dp(d(u)+(4*rho-rb)*u)-rhop*(d(u)+(4*rho-rb)*u))
             +e(ep(u)+(4*taup-tb)*u)-tau*(ep(u)+(4*taup-tb)*u)
             -(2*tau*taup-9*psi)*u)/9
        mm=mm-e(z*e(trace))/54-tau*(2*tau-tpb)*z*tilde/9
        u=z*T[n,m]
        mm=mm+2*(tau*d(u)-rho*(e(u)-4*tau*u))/9
        u=z*T[l,m]
        inner=2*e(u)+(6*tau-3*tpb)*u
        mm=mm+(2*dp(inner)-rpb*inner-(8*tau*rhop+5*rpb*tpb)*u)/9+2*e(z*M)/27

        u=z*T[l,m]
        inner=e(u)+(7*tau-4*tpb)*u
        lm=ep(inner)+(-4*taup+3*tb)*inner
        inner=d(u)+(7*rho-4*rb)*u
        lm=lm-dp(inner)-(-4*rhop+3*rpb)*inner
        lm=(lm+(-14*psi+2*psib-26*rho*rhop+32*rb*rhop-10*rpb*rb
                 -32*tau*tb+26*tau*taup+10*tb*tpb)*u)/9
        u=z*tilde
        lm=lm+(tau*d(u)-rho*e(u))/18-2*tau*(2*tau-tpb)*z*T[l,b]/9+2*rho*(2*rho-rb)*z*T[n,m]/9
        u=z*T[l,l]
        inner=2*e(u)+(6*tau-tpb)*u
        lm=lm+(2*dp(inner)+rpb*inner-(8*tau*rhop+5*rpb*tpb)*u)/18
        u=z*T[m,m]
        inner=2*d(u)+(6*rho-rb)*u
        lm=lm-(2*ep(inner)+tb*inner-(8*rho*taup+5*rb*tb)*u)/18
        lm=lm-(e(z*d(trace))+tpb*z*d(trace)+d(z*e(trace))+rb*z*e(trace))/108
        lm=lm+(d(z*M)+rb*z*M+e(z*L)+tpb*z*L)/27

        u=z*T[l,b]
        inner=e(u)+(8*tau-tpb)*u
        lb=ep(inner)-5*taup*inner
        inner=dp(u)+(4*rhop-rpb)*u
        lb=-(lb+3*(d(inner)-rho*inner)+(-6*rho*rhop+18*tau*taup)*u)/27
        lb=lb-(d(z*B)+rho*z*B+ep(z*L)+taup*z*L)/27
        u=z*T[l,l]
        inner=2*dp(u)+(9*rhop-2*rpb)*u
        lb=lb-(2*ep(inner)-taup*inner-9*rhop*taup*u)/54-2*rho*(2*rho-rb)*z*T[n,b]/9
        lb=lb+(ep(z*d(trace))+taup*z*d(trace)+d(z*ep(trace))+rho*z*ep(trace))/108
        u=z*tilde
        inner=2*d(u)+(4*rho-rb)*u
        lb=lb+(2*ep(inner)+(2*taup-3*tb)*inner-(9*rb*tb+24*rho*taup)*u)/108
        u=z*T[b,b]
        lb=lb+(rho*(e(u)+8*tau*u)-tau*(d(u)+4*rho*u))/9
        u=z*T[l,m]
        inner=ep(u)+tb*u
        lb=lb+2*(ep(inner)+(4*taup-3*tb)*inner-3*tb*taup*u)/27

        ln=-z*((2*rhop-rpb)*d(tilde)-(2*rho-rb)*dp(tilde)
               -(2*taup-tb)*e(tilde)+(2*tau-tpb)*ep(tilde))/27
        u=z*T[m,m]
        inner=ep(u)+(-tb+6*taup)*u
        ln=ln+(ep(inner)+(-2*taup+tb)*inner)/54
        u=z*T[b,b]
        inner=e(u)+(-tpb+6*tau)*u
        ln=ln-(e(inner)+(-2*tau+tpb)*inner)/54
        u=z*T[l,l]
        inner=dp(u)+(6*rhop-rpb)*u
        ln=ln-(dp(inner)+(-2*rhop+rpb)*inner)/54
        u=z*T[n,n]
        inner=d(u)+(6*rho-rb)*u
        ln=ln+(d(inner)+(-2*rho+rb)*inner)/54
        u=z*T[l,b]
        inner=dp(u)+2*rhop*u
        ln=ln-2*(e(inner)+tau*inner-(2*rpb*tpb+6*rhop*tau)*u)/27
        u=z*T[l,m]
        ln=ln-2*((2*rhop-rpb)*(ep(u)+2*taup*u)-(2*taup-tb)*(dp(u)+2*rhop*u))/27
        u=z*T[n,m]
        inner=d(u)+2*rho*u
        ln=ln+2*(ep(inner)+taup*inner-(2*rb*tb+6*rho*taup)*u)/27
        u=z*T[n,b]
        ln=ln+2*((2*rho-rb)*(e(u)+2*tau*u)-(2*tau-tpb)*(d(u)+2*rho*u))/27
        return ll,mm,lm,lb,ln
    ll,mm,lm,lb,ln=components()
    pll,pmm,plm,plb,_=components(True)
    out={(0,0):ll,(2,2):mm,(0,2):lm,(0,3):lb,(0,1):ln,(2,3):ln,
         (1,1):-pll,(3,3):-pmm,(1,3):-plm,(1,2):-plb}
    out.update({(b,a):v for (a,b),v in list(out.items())})
    dual=[[-x for x in g.cov[1]],[-x for x in g.cov[0]],g.cov[3],g.cov[2]]
    return [[sum(out[A,B].f*dual[A][a]*dual[B][b] for A in range(4) for B in range(4))
             for b in range(4)] for a in range(4)]
