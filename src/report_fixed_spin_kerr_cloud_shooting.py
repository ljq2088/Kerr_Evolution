"""Independent diagnostic complex Kerr |211> cloud by two-ended shooting.

Not imported by production modules. Angular A is analytic for complex c^2:
use general eig, never the real Hermitian eigvalsh shortcut.
"""
import argparse,hashlib,json,time
from decimal import Decimal,localcontext
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from environment_cloud import radial_coefficients
from environment_radial import horizon_series,infinity_series
ROOT=Path(__file__).resolve().parents[1]
def pair(z):return [float(z.real),float(z.imag)]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def angular_complex(ell,m,c2,size):
 ls=np.arange(abs(m),abs(m)+size+1);cosine=np.zeros((size+1,size+1))
 for j,l in enumerate(ls[:-1]):cosine[j,j+1]=cosine[j+1,j]=np.sqrt(((l+1)**2-m*m)/((2*l+1)*(2*l+3)))
 matrix=np.diag(ls[:size]*(ls[:size]+1.)).astype(complex)-complex(c2)*(cosine@cosine)[:size,:size]
 vals,vecs=np.linalg.eig(matrix);index=int(np.argmin(abs(vals-ell*(ell+1))))
 A=vals[index];v=vecs[:,index];res=float(np.linalg.norm(matrix@v-A*v)/np.linalg.norm(v))
 return A,res

def solve(a=.88,mu=.3,seed=.29629324847975716+1e-9j,outer_efolds=45.,offset=1e-5,rtol=2e-11,angular_size=20,horizon_order=4,infinity_order=6,match_radius=None):
 start=time.perf_counter();rp=1+np.sqrt(1-a*a);oh=a/(2*rp);rmin=rp+offset;rmax=2*outer_efolds/mu**2;match=2/mu**2 if match_radius is None else match_radius
 calls=0
 def mismatch(w,detail=False):
  nonlocal calls
  calls+=1;A,ares=angular_complex(1,1,a*a*(w*w-mu*mu),angular_size)
  k=1j*np.sqrt(mu*mu-w*w+0j)
  if k.imag<=0:raise ValueError('Wrong decaying branch')
  def rhs(r,y):
   d,dp,v=radial_coefficients(r,a,mu,w,1,A)
   return [y[1],-(dp*y[1]+v*y[0])/d]
  rin,din=horizon_series(rmin,a,mu,w,1,A,order=horizon_order)
  p,f,df=infinity_series(rmax,a,mu,w,1,A,k,order=infinity_order)
  options=dict(method='DOP853',rtol=rtol,atol=rtol*1e-3)
  left=solve_ivp(rhs,(rmin,match),[rin,din],**options)
  right=solve_ivp(rhs,(rmax,match),[1.+0j,1j*k+p/rmax+df/f],**options)
  if not left.success or not right.success:raise RuntimeError('Radial integration failed')
  u,v=left.y[:,-1],right.y[:,-1];res=u[1]/u[0]-v[1]/v[0]
  if detail:return dict(angular_A=pair(A),angular_eigenpair_residual=ares,log_derivative_mismatch=pair(res),left_nfev=left.nfev,right_nfev=right.nfev,left_match_state=[pair(z) for z in u],right_match_state=[pair(z) for z in v])
  return res
 def equations(x):return pair(mismatch(complex(*x)))
 initial=root(equations,pair(seed),tol=1e-11);omega=complex(*initial.x);steps=[]
 for i in range(8):
  residual=mismatch(omega)
  if abs(residual)<2e-13:break
  h=1e-6 if i<3 else 1e-7;der=(mismatch(omega+h)-mismatch(omega-h))/(2*h);correction=residual/der;steps.append(dict(mismatch=pair(residual),correction=pair(correction)));omega-=correction
 final=mismatch(omega,True)
 if abs(complex(*final['log_derivative_mismatch']))>1e-11 or not 0<omega.real<mu:raise RuntimeError(f'Invalid root: {omega}, {final}')
 return dict(a=a,mu=mu,ell=1,m=1,omega=pair(omega),Omega_H=oh,superradiant_gap=oh-omega.real,growth_or_damping='growth' if omega.imag>0 else 'damping',amplitude_timescale_M=float(1/abs(omega.imag)),m2_ionization_radius_using_real_omega=float((1/(mu-omega.real)-a)**(2/3)),shooting=dict(rmin=rmin,rmax=rmax,match_radius=match,outer_efolds=outer_efolds,horizon_offset=offset,rtol=rtol,atol=rtol*1e-3,horizon_order=horizon_order,infinity_order=infinity_order,angular_size=angular_size),initial_scipy_success=bool(initial.success),initial_scipy_message=str(initial.message),initial_nfev=int(initial.nfev),shooting_calls=calls,polish_steps=steps,final=final,elapsed_seconds=time.perf_counter()-start)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--a',type=float,default=.88);parser.add_argument('--output',type=Path,default=ROOT/'docs/environment_reproduction/fixed_spin_kerr_cloud_shooting_20260917.json');args=parser.parse_args()
 rows=[]
 for settings in [dict(outer_efolds=45.,offset=1e-5,rtol=2e-11,angular_size=20,horizon_order=4,infinity_order=6),dict(outer_efolds=70.,offset=1e-6,rtol=2e-13,angular_size=28,horizon_order=6,infinity_order=8,match_radius=25.)]:
  if rows:settings['seed']=complex(*rows[-1]['omega'])
  row=solve(a=args.a,**settings);rows.append(row);print(json.dumps(row),flush=True)
 w0,w1=(complex(*v['omega']) for v in rows)
 result=dict(status='independent_two_ended_complex_shooting_completed_not_production_cloud_replacement',method='ingoing Frobenius / decaying inverse-r boundaries + DOP853 + complex log-derivative matching; complex angular eig; no continued fraction',rows=rows,refinement_difference=dict(real=float(w1.real-w0.real),imaginary=float(w1.imag-w0.imag),absolute=float(abs(w1-w0)),imaginary_relative=float(abs((w1.imag-w0.imag)/w1.imag))),input_sha256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ROOT/'src/environment_radial.py',ROOT/'src/environment_cloud.py',ROOT/'src/environment_schwarzschild_cloud.py']},limitations=['Two numerical settings, not a certified analytic error bound.','Independent from Leaver continued fractions, but reuses local polynomial boundary expansions and primitive radial ODE coefficients.','Complex growing background is not exactly stationary or synchronized at fixed a=.88; no cloud normalization or environmental source is replaced here.'])
 refpath=ROOT/'docs/environment_reproduction/fixed_spin_cloud_frequency_20260917.json'
 if refpath.exists():
  ref=json.loads(refpath.read_text());last=ref['cases'][-1]
  if float(ref['physical_parameters']['a'])!=args.a or ref['physical_parameters']['mu'] not in ('0.3','.3'):raise ValueError('Leaver comparison has different physics')
  with localcontext() as ctx:
   ctx.prec=60
   differences=[str(Decimal(str(v))-Decimal(q)) for v,q in zip(rows[-1]['omega'],last['omega'])]
  result['independent_leaver_comparison']=dict(file=str(refpath.relative_to(ROOT)),sha256=sha(refpath),radial_terms=last['radial_terms'],omega=last['omega'],shooting_minus_leaver_reported_decimal=differences,scope='Both independent roots were computed before exchanging their numerical answers; decimal differences are not certified error bounds.')
 args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result['refinement_difference']),flush=True)

if __name__=='__main__':main()
