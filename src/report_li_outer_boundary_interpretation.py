"""Analytic interpretation of the fixed-source outer-boundary controls."""
import hashlib,json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from environment_cloud import angular_eigenvalue
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"docs/field_alignment_20260917"
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 source=OUT/"scalar22_outer_boundary_control.json";data=json.loads(source.read_text())
 li4path=OUT/"scalar22_Li_order4_boundary_control.json"
 li4=json.loads(li4path.read_text())
 control_rows=data["rows"]+li4["rows"]
 tex=ROOT/"outputs/paper_original_reference/li_2507_02045v2/source/main.tex"
 bib=tex.with_name("reference.bib")
 rows=[];physics={}
 for orbit in [41.1,42.1]:
  path=OUT/"nonstatic"/f"rp{str(orbit).replace('.','p')}_sl2_sm2_L6_q12_nr8_h32.json"
  p=json.loads(path.read_text())["parameters"];a=p["metric"]["a"];w=p["omega"];mu=p["alpha"];m=2
  k2=w*w-mu*mu;beta=2*w*w-mu*mu;A=angular_eigenvalue(2,2,a*a*k2);C=A+a*a*w*w-2*a*m*w
  C2=12*w*w-4*mu*mu-a*a*k2-A
  def approximate(r):return k2+2*beta/r+C2/r**2
  def exact(r):
   delta=r*r-2*r+a*a;K=(r*r+a*a)*w-a*m
   return (K*K+1-a*a)/delta**2-(mu*mu*r*r+C)/delta
  roots=np.roots([k2,2*beta,C2])
  row=dict(orbit=orbit,a=a,mu=mu,omega=w,k_squared=k2,beta=beta,
   angular_A=float(A),C2=float(C2),long_range_scale_2beta_over_abs_k2=2*beta/abs(k2),
   approximate_turning_points=sorted(float(r.real) for r in roots if abs(r.imag)<1e-12 and r.real>0),
   tests=[dict(radius=r,Q_truncated=float(approximate(r)),Q_exact=float(exact(r)),
     region="oscillatory" if exact(r)>0 else "evanescent",
     beta_over_abs_k2_r=beta/(abs(k2)*r)) for r in [1000.,4000.,8000.,32000.,64000.]])
  if k2<0:row["exact_outer_turning_point"]=brentq(exact,4000.,8000.,xtol=1e-8)
  else:
   row["leading_order_current_normalization"]=[]
   for r in [1000.,4000.,8000.]:
    observed=next(x for x in data["rows"] if x["orbit"]==orbit and x["inverse_r_order"]==0 and x["outer_radius"]==r)["actual_current_over_nominal_FI"]
    predicted=(1-2/r+a*a/r**2)*(1+beta/(k2*r))
    row["leading_order_current_normalization"].append(dict(radius=r,analytic_factor=predicted,numerical_factor=observed,relative_difference=observed/predicted-1))
  rows.append(row);physics[orbit]=(k2,beta,C2)
 result=dict(status="analytic_interpretation_of_deliberately_finite_boundary_controls",
   transform="y=sqrt(Delta) R; y''+Q y=0",
   exact_Q="Q=(K^2+1-a^2)/Delta^2-(mu^2 r^2+A+a^2 omega^2-2 a m omega)/Delta",
   expansion="Q=k^2+2 beta/r+C2/r^2+O(r^-3), beta=2 omega^2-mu^2, C2=12 omega^2-4 mu^2-a^2 k^2-A",
   leading_current_ratio="N_actual/(2 k |ZI|^2)=(Delta/Rout^2)*(1+beta/(k^2 Rout)), for R_up=exp(i k r) r^(-1+i beta/k) exp(-2ik log2)",
   finite_order_warning="The 0/6 controls use the leading asymptotic prefactor. A separate control correctly uses the published exact tortoise prefactor with fourth-order coefficients independently derived from the ODE.",
   rows=rows,source_sha256=sha(source),li_exact_order4_result_sha256=sha(li4path),
   Li_method_evidence=dict(main_tex=str(tex.relative_to(ROOT)),main_tex_sha256=sha(tex),
     boundary_equation_lines=[568,575],expansion_lines=[592,604],specified_highest_order=4,
     cutoff_comment_line=729,explicit_outer_cutoff_found_in_local_main_tex=False,
     supplementary_boundary_repository="https://github.com/dongjun826/EMRI-in-scalar-clouds.git",
     bibliography_file=str(bib.relative_to(ROOT)),bibliography_sha256=sha(bib),
     conclusion="Expansion order is known. This local text search does not establish the numerical outer endpoint or the actual implementation parameters."),
   interpretation="Finite external asymptotics can alter the bound-mode nodes and the formal radiation normalization. This demonstrates a plausible failure class, not the original authors' actual boundary choice or a proved explanation of their plotted discrepancy.",
   implementation_sha256=sha(__file__))
 (OUT/"scalar22_outer_boundary_interpretation.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
 fig,axes=plt.subplots(2,2,figsize=(12,8),layout="constrained")
 r=np.geomspace(500,15000,400);k2,beta,c2=physics[42.1]
 ax=axes[0,0];ax.plot(r,(k2+2*beta/r+c2/r**2)*1e5,color="black")
 for radius in [4000,8000]:ax.axvline(radius,ls=":",label=f"Boundary {radius}M")
 ax.axhline(0,color="0.6",lw=.8);ax.set(xscale="log",xlabel="r/M",ylabel="Q x 1e5",title="Bound branch: outer turning point ~5170M");ax.legend()
 ax=axes[0,1]
 for order,color in [(0,"C0"),(4,"C2"),(6,"C1")]:
  sub=[x for x in control_rows if x["orbit"]==41.1 and x["inverse_r_order"]==order]
  x=[v["outer_radius"] for v in sub]
  ax.plot(x,[v["total_infinity_flux_only_scalar22_replaced"]["nominal_total_over_Li"] for v in sub],"o-",color=color,label=f"Order {order}: nominal infinity")
  ax.plot(x,[v["total_infinity_flux_only_scalar22_replaced"]["actual_current_total_over_Li"] for v in sub],"s--",color=color,label=f"Order {order}: conserved current")
 ax.axhline(1.2182028725,color="black",ls=":",label="Coulomb 32k");ax.axhline(1,color=".6",lw=.8)
 ax.set(xscale="log",xlabel="Outer boundary / M",ylabel="Total flux / Li vector reference",title="Orbit 41.1M: nominal and actual flux differ");ax.legend(fontsize=7)
 for i,radius in enumerate([50,150]):
  ax=axes[1,i]
  for order,color in [(0,"C0"),(4,"C2"),(6,"C1")]:
   sub=[x for x in control_rows if x["orbit"]==42.1 and x["inverse_r_order"]==order]
   errors=[next(t["normalized_amplitude_rms_difference"] for t in v["full18_only_scalar22_replaced"] if t["radius"]==radius) for v in sub]
   ax.plot([v["outer_radius"] for v in sub],errors,"o-",color=color,label=f"Order {order}")
  base=next(x for x in data["rows"] if x["orbit"]==42.1 and x["boundary_method"]=="coulomb" and x["outer_radius"]==32000)
  error=next(t["normalized_amplitude_rms_difference"] for t in base["full18_only_scalar22_replaced"] if t["radius"]==radius)
  ax.axhline(error,color="black",ls=":",label="Coulomb 32k")
  ax.set(xscale="log",yscale="log",xlabel="Outer boundary / M",ylabel="Angular amplitude RMS difference / reference RMS",title=f"Orbit 42.1M, r={radius}M: only scalar22 changed")
  ax.legend()
 fig.suptitle("Frozen sources; finite-boundary sensitivity is not author-code identification",fontsize=12)
 for suffix in ["png","pdf"]:fig.savefig(OUT/f"scalar22_outer_boundary_control.{suffix}",dpi=160)
 plt.close(fig)
 print(json.dumps(rows,indent=2))
if __name__=="__main__":main()
