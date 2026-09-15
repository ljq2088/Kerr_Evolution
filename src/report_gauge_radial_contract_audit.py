"""Check spin 0/1 radial conventions and Taylor recurrence independently.

This does not change production. It tests exact symmetry/TS contracts rather
than comparing HBL and TEUK, which share boundary and spin-flip code.
"""
from pathlib import Path
import hashlib, json, math
import importlib
import mpmath as mp
import numpy as np
from pybhpt.radial import RadialTeukolsky
from lorenz_ghp import KerrGHP
from lorenz_jet import Jet
from lorenz_mode_jet import separated_jet
from environment_cloud import angular_eigenvalue

ROOT=Path(__file__).resolve().parents[1]
A=.8771530275949366
R0=20.
RP=1+np.sqrt(1-A*A)
RADII=np.array([RP+.0005,2.,10.,19.9,20.1,40.,120.,316.])


def encode(z):
    a=np.asarray(z,dtype=complex)
    return np.stack([a.real,a.imag],axis=-1).tolist()


def solve(s,l,m):
    omega=m/(R0**1.5+A)
    q=RadialTeukolsky(s,l,m,A,omega,RADII)
    q.solve()
    return q,np.array([q(b,d) for b in ('In','Up') for d in (0,1)]).reshape(2,2,-1)


def polynomial_taylor(s,lam,m,omega,r,value,derivative,n=8):
    """80-digit recurrence from polynomial Delta^2 R''+B R'+C R=0.

    Does not call Jet or differentiate the divided ODE used in production.
    The double-precision input state is held fixed, so this isolates recurrence
    roundoff and does not establish accuracy of the initial radial state.
    """
    with mp.workdps(80):
        aa=mp.mpf(A);rr=mp.mpf(r);ww=mp.mpf(omega);ll=mp.mpf(lam)
        delta=lambda x:x*x-2*x+aa*aa
        kk=lambda x:(x*x+aa*aa)*ww-aa*m
        functions=(lambda x:delta(x)**2,
            lambda x:(s+1)*delta(x)*2*(x-1),
            lambda x:kk(x)**2-2j*s*(x-1)*kk(x)+delta(x)*(4j*s*ww*x-ll))
        poly=[mp.taylor(f,rr,4) for f in functions]
        c=[mp.mpc(value),mp.mpc(derivative)]
        for k in range(n-1):
            total=mp.mpc(0)
            for degree,v in enumerate(poly[0][1:],1):
                j=k-degree+2
                if 2<=j<len(c):total+=v*j*(j-1)*c[j]
            for degree,v in enumerate(poly[1]):
                j=k-degree+1
                if 1<=j<len(c):total+=v*j*c[j]
            for degree,v in enumerate(poly[2]):
                j=k-degree
                if 0<=j<len(c):total+=v*c[j]
            c.append(-total/(poly[0][0]*(k+2)*(k+1)))
        return np.array([complex(z) for z in c])


def source_normalization_audit():
    """Independent complex changes of In/Up basis must cancel in each field.

    Includes both Maxwell chiralities, their sourced amplitudes and the compact
    scalar chi. Trace and kappa use a different solver and are held fixed.
    """
    import lorenz_metric
    modules=[importlib.import_module(name) for name in
             ('lorenz_metric','lorenz_spin1','lorenz_spin1_chiral','lorenz_chi')]
    def clear():
        for module in modules:
            for f in vars(module).values():
                if callable(getattr(f,'cache_clear',None)):f.cache_clear()
    def values():
        out={}
        for r in (10.,40.):
            for kind in (0,1):
                if kind==1:
                    _,h=lorenz_metric.spin1_metric(r,1.1,R0,A,1,1,order=6,full_current=True)
                else:
                    _,h=lorenz_metric.spin0_metric(r,1.1,R0,A,1,1,order=6)
                out[(r,kind)]=np.array([[v.value for v in row] for row in h],complex)
        return out
    clear();baseline=values()
    class ScaledRadial(RadialTeukolsky):
        def factor(self,bc):
            sign=1 if self.azimuthalmode>0 else -1
            return ((.7+.13*self.spinweight)+1j*(.4+.11*sign) if bc=='In'
                    else (-1.3+.17*self.spinweight)+1j*(.2-.07*sign))
        def radialsolution(self,bc,pos):
            return self.factor(bc)*super().radialsolution(bc,pos)
        def radialderivative(self,bc,pos):
            return self.factor(bc)*super().radialderivative(bc,pos)
        def radialderivative2(self,bc,pos):
            return self.factor(bc)*super().radialderivative2(bc,pos)
    originals=[module.RadialTeukolsky for module in modules]
    try:
        for module in modules:module.RadialTeukolsky=ScaledRadial
        clear();changed=values()
        rows=[]
        for key,b in baseline.items():
            c=changed[key]
            rows.append(dict(r=key[0],sector='full_spin1_both_chiralities' if key[1] else 'spin0_trace_kappa_chi',
                relative_metric_frobenius_change=float(np.linalg.norm(c-b)/np.linalg.norm(b)),
                max_absolute_component_change=float(np.max(abs(c-b)))))
        return rows
    finally:
        for module,original in zip(modules,originals):module.RadialTeukolsky=original
        clear()


def main():
    conjugation=[];ts=[];jets=[];eigenvalues=[]
    for ell,m in [(1,1),(6,1),(18,1),(6,6)]:
        omega=m/(R0**1.5+A)
        positive={};negative={}
        for s in (0,-1,1):
            q,v=solve(s,ell,m);qn,vn=solve(s,ell,-m)
            positive[s]=(q,v);negative[s]=(qn,vn)
            err=abs(vn/v.conjugate()-1)
            physical=np.where(RADII[None,:]<R0,err[0],err[1])
            conjugation.append(dict(spin=s,ell=ell,m=m,lambda_difference=float(qn.eigenvalue-q.eigenvalue),
                max_relative_all_branches=float(np.max(err)),
                physical_max_relative=float(np.max(physical)),
                relative_by_branch_derivative_radius=err.tolist()))
            if s==0:
                scalar_A=angular_eigenvalue(ell,m,A*A*omega*omega)
                eigenvalues.append(dict(ell=ell,m=m,pybhpt_radial_lambda=q.eigenvalue,
                  scalar_angular_A=scalar_A,
                  lambda_shift_residual=float(q.eigenvalue-(scalar_A+A*A*omega*omega-2*A*m*omega))))
            for index in (1,3,4,7):
                r=RADII[index];side=0 if r<R0 else 1
                v0,v1=v[side,:,index]
                reference=polynomial_taylor(s,q.eigenvalue,m,omega,float(r),v0,v1)
                for dtype in (np.complex128,np.clongdouble):
                    Jet.coefficient_dtype=dtype
                    g=KerrGHP(float(r),1.1,A,omega=omega,m=m,order=8)
                    # Setting S0=1 and S1=0 isolates the radial slice. Undo the
                    # Kinnersley zeta factor before extracting Taylor values.
                    field=separated_jet(g,s,q.eigenvalue,v0,v1,1.,0.)
                    raw=field*g.zeta**(abs(s)-s)
                    got=np.array(raw.c[:,0],dtype=complex)
                    relative=abs((got-reference)/reference)
                    jets.append(dict(spin=s,ell=ell,m=m,r=float(r),bc=('In','Up')[side],
                      coefficient_dtype=np.dtype(dtype).name,
                      max_relative_derivative_orders_2_to_8=float(np.max(relative[2:])),
                      relative_by_derivative_order=relative.tolist()))
        Jet.coefficient_dtype=np.complex128
        minus=positive[-1];plus=positive[1]
        lam=minus[0].eigenvalue
        delta=RADII*RADII-2*RADII+A*A;dp=2*(RADII-1)
        K=(RADII*RADII+A*A)*omega-A*m
        q=K/delta;qp=2*RADII*omega/delta-K*dp/delta**2
        Kh=2*RP*omega-A*m;hgap=2*np.sqrt(1-A*A)
        B2=lam*lam+4*A*m*omega-4*A*A*omega*omega
        for side,bc in enumerate(('In','Up')):
            R,Rp=minus[1][side]
            V=(K*K+2j*(RADII-1)*K)/delta-4j*omega*RADII-lam
            Rpp=-V*R/delta
            transformed=Rpp-2j*q*Rp-(1j*qp+q*q)*R
            constant=-4*Kh*Kh-2j*Kh*hgap if bc=='In' else -B2/(4*omega*omega)
            error=abs(transformed/(constant*plus[1][side,0])-1)
            used=RADII<R0 if bc=='In' else RADII>R0
            ts.append(dict(ell=ell,m=m,bc=bc,lambda_plus_minus_residual=float(plus[0].eigenvalue-lam+2),
                predicted_TS_constant=encode(constant),
                measured_TS_ratio=encode(transformed/plus[1][side,0]),
                max_relative_all_radii=float(np.max(error)),physical_max_relative=float(np.max(error[used])),
                relative_by_radius=error.tolist()))
    normalization=source_normalization_audit()
    paths=[Path(__file__),ROOT/'src/lorenz_metric.py',ROOT/'src/lorenz_mode_jet.py',ROOT/'src/lorenz_spin1.py',
      ROOT/'src/lorenz_spin1_chiral.py',ROOT/'src/lorenz_chi.py',ROOT/'src/lorenz_kappa.py',ROOT/'src/environment_radial.py']
    output=dict(status='radial_convention_and_recurrence_audit_not_independent_flux_validation',a=A,rp=R0,radii=RADII.tolist(),
      source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
      radial_equation="Delta R''+(s+1)Delta'R'+[(K^2-2is(r-1)K)/Delta+4is omega r-lambda]R=0",
      radial_lambda='A_s + a^2 omega^2 - 2 a m omega; lambda_+1=lambda_-1-2',
      conjugation_identity='R_s(-m,-omega)=conj(R_s(m,omega)) for the same retarded In/Up normalization',
      ts_identity="(d/dr-iK/Delta)^2 R_-1 = C_bc R_+1",
      ts_normalization="C_In=-4 K_+^2-2i K_+(r_+-r_-); C_Up=-B^2/(4 omega^2); B^2=lambda_-1^2+4amomega-4a^2omega^2",
      conjugation=conjugation,eigenvalues=eigenvalues,spin1_ts=ts,taylor_recurrence=jets,
      source_basis_normalization_invariance=normalization,
      limitations=['TS may share algebraic identities with the library spin flip; boundary constants are derived here, not fitted.',
        'Conjugation can pass with a shared real normalization defect.',
        'Polynomial Taylor reference holds double radial values fixed; it tests the recurrence, not radial input accuracy.',
        'Trace/kappa RadialGreen paths are not changed by the pybhpt control and require separate mass-derivative and boundary checks.'])
    dest=ROOT/'docs/environment_reproduction/gauge_radial_contract_audit.json'
    dest.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(dict(conjugation_max=max(z['physical_max_relative'] for z in conjugation),
      ts_max=max(z['physical_max_relative'] for z in ts),
      lambda_max=max(abs(z['lambda_shift_residual']) for z in eigenvalues),
      source_basis_normalization_invariance_max=max(z['relative_metric_frobenius_change'] for z in normalization),
      taylor_double_max=max(z['max_relative_derivative_orders_2_to_8'] for z in jets if z['coefficient_dtype']=='complex128'),
      taylor_extended_max=max(z['max_relative_derivative_orders_2_to_8'] for z in jets if z['coefficient_dtype']!='complex128'))))

if __name__=='__main__':main()
