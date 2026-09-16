"""Li exact tortoise-prefactor, fourth-order boundary: frozen-source control."""
import hashlib,json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import mpmath as mp
import numpy as np
import report_li_scalar22_outer_boundary_control as base
ROOT=base.ROOT;OUT=base.OUT;RADII=base.RADII
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def formal_boundary(a,mu,w,m,A,r,order=4,dps=65,propagating=True):
 with mp.workdps(dps):
  a,mu,w,A,r=[mp.mpf(str(float(v))) for v in (a,mu,w,A,r)]
  k=mp.sqrt(w*w-mu*mu) if propagating else 1j*mp.sqrt(mu*mu-w*w)
  size=order+4
  def series(values):
   out=[mp.mpc(0)]*size
   for i,v in enumerate(values):out[i]=mp.mpc(v)
   return out
  def add(*items):return [sum(x[i] for x in items) for i in range(size)]
  def scale(x,c):return [v*c for v in x]
  def mul(x,y):return [sum(x[j]*y[i-j] for j in range(i+1)) for i in range(size)]
  def inverse(x):
   y=[1/x[0]]
   for i in range(1,size):y.append(-sum(x[j]*y[i-j] for j in range(1,i+1))/x[0])
   return y
  def shift(x,n):return [mp.mpc(0)]*n+x[:size-n]
  f=series([1,-2,a*a]);g=series([1,0,a*a])
  # Exact logarithmic derivative of the paper's P(r).
  L=add(scale(mul(g,inverse(f)),1j*k),series([0,1j*mu*mu/k]),scale(shift(inverse(g),1),-1))
  Lprime=series([0,0]+[-(i-1)*L[i-1] for i in range(2,size)])
  aa=mul(f,f);bb=series([0,2,-6,4+2*a*a,-2*a*a])
  b=a*a*w-a*m;C=A+a*a*w*w-2*a*m*w
  cc=series([k*k,2*mu*mu,2*w*b-C-a*a*mu*mu,2*C,b*b-a*a*C])
  E=add(mul(aa,add(Lprime,mul(L,L))),mul(bb,L),cc)
  D=add(scale(shift(mul(aa,L),2),-2),scale(shift(bb,2),-1),scale(shift(aa,3),2))
  def coefficient(n,power):
   def at(arr,i):return arr[i] if 0<=i<size else 0
   return at(E,power-n)+n*at(D,power-n+1)+n*(n-1)*at(aa,power-n-2)
  if max(abs(E[0]),abs(E[1]))>mp.mpf("1e-50"):raise ValueError("Leading asymptotic cancellation failed")
  coeff=[mp.mpc(1)]
  for n in range(1,order+1):
   power=n+1
   coeff.append(-sum(coeff[j]*coefficient(j,power) for j in range(n))/coefficient(n,power))
  residual=[sum(coeff[j]*coefficient(j,power) for j in range(len(coeff))) for power in range(order+2)]
  B=sum(v/r**n for n,v in enumerate(coeff))
  Bd=sum(-n*v/r**(n+1) for n,v in enumerate(coeff) if n)
  rp=1+mp.sqrt(1-a*a);rm=2-rp;delta=(r-rp)*(r-rm)
  rstar=r+2*rp/(rp-rm)*mp.log((r-rp)/2)-2*rm/(rp-rm)*mp.log((r-rm)/2)
  P=mp.exp(1j*k*rstar)*mp.exp(1j*mu*mu/k*mp.log(r))/mp.sqrt(r*r+a*a)
  logP=1j*k*(r*r+a*a)/delta+1j*mu*mu/(k*r)-r/(r*r+a*a)
  up=P*B if propagating else mp.mpc(1)
  derivative=up*(logP+Bd/B)
  encode=lambda z:[mp.nstr(mp.re(z),45),mp.nstr(mp.im(z),45)]
  info=dict(order=order,precision_digits=dps,coefficients=[encode(z) for z in coeff],
    recurrence_residual_through_power_5=[encode(z) for z in residual],
    relative_last_term=float(abs(coeff[-1]/r**order)/abs(B)),
    partial_sum_B=encode(B),boundary_log_derivative=encode(logP+Bd/B),
    exact_prefactor="exp(i k r*) r^(i mu^2/k)/sqrt(r^2+a^2)",
    construction="Independent formal series from P'/P and polynomial radial ODE, no use of production infinity_series coefficients.")
  return complex(up),complex(derivative),info

def calculate(orbit_and_outer):
 import environment_radial as radial
 orbit,outer=orbit_and_outer
 original=radial.coulomb_boundary;info={}
 def paper_boundary(r,a,mu,omega,lam,k,propagating):
  u,du,audit=formal_boundary(a,mu.real if hasattr(mu,"real") else mu,omega,2,lam,r,propagating=propagating)
  info.update(audit)
  return u,du
 try:
  radial.coulomb_boundary=paper_boundary
  row=base.calculate((orbit,"coulomb",None,outer))
 finally:radial.coulomb_boundary=original
 row.update(boundary_method="Li_exact_prefactor_series4",inverse_r_order=4,
   boundary_status="preselected_finite_Rout_control_not_author_endpoint_identification",
   paper_boundary_audit=info)
 return row

def main():
 from environment_wake import EnvironmentalWake
 from source_provenance import source_fingerprint,validate_saved_samples
 current=source_fingerprint()
 oldpath=OUT/"scalar22_outer_boundary_control.json";old=json.loads(oldpath.read_text())
 with ProcessPoolExecutor(4) as pool:
  rows=list(pool.map(calculate,[(orbit,r) for orbit in [41.1,42.1] for r in [4000.,8000.,32000.]]))
 # Save the new independent coefficients early for the author-notebook audit.
 out=OUT/"scalar22_Li_order4_boundary_control.json"
 preliminary=dict(status="radial_computations_complete_field_comparison_pending",rows=rows)
 out.write_text(json.dumps(preliminary,indent=2,allow_nan=False)+"\n")
 alignmentpath=OUT/"li_field_alignment.json";alignment=json.loads(alignmentpath.read_text())
 fluxpath=OUT/"li_new_radius_flux.json";fluxreference=json.loads(fluxpath.read_text())
 polarpath=ROOT/"docs/environment_reproduction/li_field_polar_reference_20260917.npz";polar=np.load(polarpath)
 inputs={str(p.relative_to(ROOT)):sha(p) for p in [oldpath,alignmentpath,fluxpath,polarpath,Path(__file__),ROOT/"src/report_li_scalar22_outer_boundary_control.py"]}
 for orbit in [41.1,42.1]:
  reports=[]
  for relative,h in alignment["inputs_sha256"].items():
   p=ROOT/relative
   if p.suffix!=".json":continue
   d=json.loads(p.read_text());param=d.get("parameters",{})
   if param.get("metric",{}).get("orbital_radius")!=orbit or "samples" not in d:continue
   assert sha(p)==h;validate_saved_samples(d,current);reports.append(d);inputs[relative]=h
  assert len(reports)==18
  wake=EnvironmentalWake(reports,ellmax=5,allow_mixed_discretization=True)
  field=wake.evaluate(RADII[:,None],np.pi/2,polar["phi"][None,:])/.3**3
  baseline=next(r for r in old["rows"] if r["orbit"]==orbit and r["boundary_method"]=="coulomb" and r["outer_radius"]==32000)
  ref=next(r for r in fluxreference["rows"] if r["orbit"]==orbit)
  old22=next(m["orbital_energy_infinity"] for m in ref["modes"] if (m["ell"],m["m"])==(2,2))
  othertotal=ref["local_infinity_flux"]-old22;Li=ref["references"]["Li"]["reference"]
  refindex=int(np.flatnonzero(polar["rp"]==orbit)[0])
  for row in [r for r in rows if r["orbit"]==orbit]:
   z=base.dec(row["scalar22_equatorial_per_q_epsilon"]);zb=base.dec(baseline["scalar22_equatorial_per_q_epsilon"])
   row["scalar22_amplitude_ratio_to_coulomb32k"]=(abs(z)/abs(zb)).tolist()
   row["scalar22_complex_relative_difference_to_coulomb32k"]=(abs(z-zb)/abs(zb)).tolist()
   changed=field+(z-zb)[:,None]*np.exp(2j*polar["phi"])[None,:]
   row["full18_only_scalar22_replaced"]=[]
   for i,radius in enumerate(RADII):
    inds=np.flatnonzero(np.isclose(polar["radii"],radius,rtol=0,atol=1e-10))
    amp=abs(changed[i])
    if not len(inds):
     row["full18_only_scalar22_replaced"].append(dict(radius=radius,reference_available=False,local_angular_rms=float(np.sqrt(np.mean(amp**2)))));continue
    j=int(inds[0]);values=polar["nominal_abs_field"][refindex,j]
    lo=polar["lower_color_interval"][refindex,j];hi=polar["upper_color_interval"][refindex,j]
    rms=float(np.sqrt(np.mean(values**2)))
    row["full18_only_scalar22_replaced"].append(dict(radius=radius,local_angular_rms=float(np.sqrt(np.mean(amp**2))),reference_nominal_angular_rms=rms,angular_rms_ratio=float(np.sqrt(np.mean(amp**2))/rms),normalized_amplitude_rms_difference=float(np.sqrt(np.mean((amp-values)**2))/rms),fraction_inside_raster_color_interval=float(np.mean((amp>=lo)&(amp<=hi)))))
   nominal=othertotal+row["nominal_asymptotic_flux"]["infinity"]["orbital_energy"]
   actual=othertotal+row["source_free_actual_orbital_current"][0]
   row["total_infinity_flux_only_scalar22_replaced"]=dict(unchanged_other_channels=othertotal,nominal_asymptotic_total=nominal,actual_current_total=actual,Li_vector_plot_reference=Li,nominal_total_over_Li=nominal/Li,actual_current_total_over_Li=actual/Li)
 if source_fingerprint()!=current:raise ValueError("Source fingerprint changed")
 result=dict(status="paper_prefactor_fourth_order_fixed_source_control_completed",
  rows=rows,inputs_sha256=inputs,source_fixed=True,production_files_modified=False,
  source_provenance=current,
  scope="Correct published prefactor and fourth-order formal recurrence; independent coefficient derivation. Numerical outer endpoints are preselected controls, not known author values.",
  limitations=["Only scalar22 propagator changed; all sources and 17 other modes frozen.",
   "Author notebook coefficients require their own independent comparison before claiming original-author code execution.",
   "Near-threshold last-term/endpoint convergence must be checked; finite Wronskian residual is insufficient.",
   "Local source still has stationary a=.877153; this is not a newly computed fixed-a=.88 cloud/source."])
 out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
 for r in rows:print(r["orbit"],r["outer_radius"],"lastterm",r["paper_boundary_audit"]["relative_last_term"],"A ratios",r["scalar22_amplitude_ratio_to_coulomb32k"],"nominal/actual total vs Li",r["total_infinity_flux_only_scalar22_replaced"],"ring errors",[v.get("normalized_amplitude_rms_difference") for v in r["full18_only_scalar22_replaced"]],flush=True)
 print(out)
if __name__=="__main__":main()
