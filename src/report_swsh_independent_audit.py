"""Cross-check pybhpt against local qnm's independently assembled angular matrix.

qnm decorators alone are removed to run its Python formulas without numba.
Harmonic values/derivatives are additionally checked using 50-digit Wigner sums.
"""
import ast, hashlib, json, math
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.linalg import eigh
from pybhpt.swsh import SpinWeightedSpheroidalHarmonic
from pybhpt.radial import RadialTeukolsky


def load_reference(path):
    raw=path.read_bytes();tree=ast.parse(raw)
    tree.body=[n for n in tree.body if not(isinstance(n,ast.ImportFrom) and n.module=='numba')]
    for n in ast.walk(tree):
        if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef)):n.decorator_list=[]
    ns={};exec(compile(tree,str(path),'exec'),ns)
    return ns,hashlib.sha256(raw).hexdigest()


def spherical(s,l,m,t):
    pref=(-1)**s*mp.sqrt(mp.mpf(2*l+1)/(4*mp.pi)*
        math.factorial(l-s)*math.factorial(l+s)*math.factorial(l-m)*math.factorial(l+m))
    return pref*sum((-1)**(m+s+k)*mp.cos(t/2)**(2*l-s-m-2*k)*mp.sin(t/2)**(m+s+2*k)/
        (math.factorial(l-s-k)*math.factorial(k)*math.factorial(m+s+k)*math.factorial(l-m-k))
        for k in range(max(0,-m-s),min(l-s,l-m)+1))


def main():
    reference=Path('/home/ljq/code/qnm/qnm/angular.py');q,sha=load_reference(reference)
    mp.mp.dps=50;rows=[];a=.8771530275949366
    for r0 in (10.,20.,41.8):
      for m in (-1,1):
       omega=m/(r0**1.5+a);c=a*omega
       for s in (-2,-1,0,1,2):
        low=max(abs(s),abs(m))
        for ell in sorted({low,6,18}):
          h=SpinWeightedSpheroidalHarmonic(s,ell,m,c);n=len(h.coeffs)
          matrix=np.array([[q['M_matrix_elem'](s,c,m,l,k) for k in range(low,low+n)] for l in range(low,low+n)])
          vals,vec=eigh(matrix);idx=ell-low;coeff=vec[:,idx]*np.sign(vec[idx,idx])
          # qnm returns A; pybhpt uses lambda=A+c^2-2mc.
          lam=vals[idx]+c*c-2*m*c
          radial=RadialTeukolsky(s,ell,m,a,omega,np.array([r0+1.]))
          row=dict(r0=r0,s=s,ell=ell,m=m,c=c,
            matrix_max_error=float(np.max(abs(h.sparse_matrix(n).toarray()-matrix-(c*c-2*m*c)*np.eye(n)))),
            lambda_absolute_error=float(abs(h.eigenvalue-lam)),
            radial_lambda_absolute_error=float(abs(radial.eigenvalue-lam)),
            coefficient_l2_error=float(np.linalg.norm(h.coeffs-coeff)))
          if r0==20. and ell==low:
            theta=[.1,.7,1.1,np.pi/2,2.7,3.04];errors=[]
            # Independent Wigner evaluation retains every reference coefficient.
            def f(t):return sum(mp.mpf(float(b))*spherical(s,low+j,m,t) for j,b in enumerate(coeff))
            for derivative in (0,1,2):
                expected=np.array([float(mp.diff(f,mp.mpf(float(t)),derivative)) for t in theta])
                actual=h(np.array(theta),deriv=derivative)
                errors.append(float(max(abs(actual-expected))/max(abs(expected))))
            row['wigner_relative_max_errors_value_dtheta_dtheta2']=errors
          rows.append(row)
      print('Completed orbit',r0,flush=True)
    from environment_source import angular_mode
    scalar_rows=[]
    wc=.29629324847975713;mu=.3;r0=20.
    for ell,m in ((1,1),(0,0),(2,2),(5,-5),(6,6)):
        omega=wc+(m-1)/(r0**1.5+a);c2=a*a*(omega*omega-mu*mu)
        c=complex(np.lib.scimath.sqrt(c2));low=abs(m);n=32
        matrix=np.array([[q['M_matrix_elem'](0,c,m,l,k) for k in range(low,low+n)] for l in range(low,low+n)])
        assert np.max(abs(matrix.imag))<1e-14
        vals,vec=eigh(matrix.real);idx=ell-low;coef=vec[:,idx]*np.sign(vec[idx,idx])
        theta=np.array([.1,.7,1.1,np.pi/2,2.7,3.04])
        value,derivative,lam=angular_mode(theta,ell,m,c2)
        def f(t):return sum(mp.mpf(float(b))*spherical(0,low+j,m,t) for j,b in enumerate(coef))
        expected=np.array([float(f(mp.mpf(float(t)))) for t in theta])
        dexpected=np.array([float(mp.diff(f,mp.mpf(float(t)))) for t in theta])
        x,w=np.polynomial.legendre.leggauss(96)
        scalar_rows.append(dict(ell=ell,m=m,c_squared=c2,
            A_absolute_error=float(abs(lam-vals[idx])),
            value_relative_max_error=float(max(abs(value-expected))/max(abs(expected))),
            derivative_relative_max_error=float(max(abs(derivative-dexpected))/max(abs(dexpected))),
            full_sphere_norm=float(2*np.pi*np.dot(w,angular_mode(np.arccos(x),ell,m,c2)[0]**2))))
    result=dict(status='independent_angular_cross_check_not_full_flux_validation',scalar_rows=scalar_rows,reference=str(reference),
      reference_sha256=sha,reference_execution='AST removes numba import and decorators only',
      convention='lambda_pybhpt = A_qnm + (a omega)^2 - 2m a omega',rows=rows)
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/swsh_independent_qnm_audit.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:max(row[k] for row in rows) for k in ('matrix_max_error','lambda_absolute_error','radial_lambda_absolute_error','coefficient_l2_error')},indent=2))
    print('Wigner maxima',np.max([r['wigner_relative_max_errors_value_dtheta_dtheta2'] for r in rows if 'wigner_relative_max_errors_value_dtheta_dtheta2' in r],axis=0))


if __name__=='__main__':main()
