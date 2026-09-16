"""Frozen fresh scalar22 source: deliberately finite asymptotic boundaries."""
import hashlib,json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"docs/field_alignment_20260917"
RADII=np.array([20.,50.,100.,150.,200.])
CONFIGS=[("coulomb",None,32000.),("coulomb",None,64000.)]+[("series",n,r) for n in [0,6] for r in [1000.,4000.,8000.]]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def enc(z):
 z=np.asarray(z,complex);return np.stack([z.real,z.imag],axis=-1).tolist()
def dec(z):
 z=np.asarray(z);return z[...,0]+1j*z[...,1]
def calculate(job):
 from environment_cloud import mode_flux
 import environment_radial as radial
 from environment_response import SampledResponse
 from environment_source import angular_mode
 from source_provenance import source_fingerprint,validate_saved_samples
 orbit,method,order,outer=job
 path=OUT/"nonstatic"/f"rp{str(orbit).replace('.','p')}_sl2_sm2_L6_q12_nr8_h32.json"
 data=json.loads(path.read_text());validate_saved_samples(data,source_fingerprint())
 p=data["parameters"];a=p["metric"]["a"];mu=p["alpha"];w=p["omega"];wc=w-1/(orbit**1.5+a)
 original=radial.infinity_series
 # A local override in one worker process; never changes a production file.
 if method=="series":
  def fixed_order(r,a,mu,omega,m,lam,k,order=6):
   return original(r,a,mu,omega,m,lam,k,order=forced_order)
  forced_order=order
  radial.infinity_series=fixed_order
 started=time.perf_counter()
 try:
  green=radial.RadialGreen(a,mu,w,2,2,rmax=outer,offset=p["green_horizon_offset"],rtol=1e-11,infinity_method=method)
 finally:
  radial.infinity_series=original
 samples=data["samples"];r=np.array([x["r"] for x in samples]);J=np.array([complex(*x["source"]) for x in samples]);weights=np.array([x["weight"] for x in samples])
 response=SampledResponse(green,p["source_panels"],r,J,log_first=True)
 R,dR=response.evaluate(RADII)
 zi=response.up_coefficient;zh=response.horizon_coefficient
 S=angular_mode(np.pi/2,2,2,a*a*(w*w-mu*mu))[0]
 flux=mode_flux(w,2,wc,1,mu,a,zi,zh)
 checkr=np.array([400.,min(800.,outer)])
 field=response.evaluate(checkr);delta=checkr**2-2*checkr+a*a
 actualcharge=2*delta*np.imag(np.conjugate(field[0])*field[1])
 actualenergy=(w-wc)*actualcharge
 _,f6,_=original(outer,a,mu,w,2,green.lam,green.k,order=6)
 _,f5,_=original(outer,a,mu,w,2,green.lam,green.k,order=5)
 gausszi=np.dot(weights,green.insol.sol(r)[0]*J)/green.w0
 gausszh=np.dot(weights,green.upsol.sol(r)[0]*J)/green.w0
 return dict(orbit=orbit,boundary_method=method,inverse_r_order=order,outer_radius=outer,
    boundary_status=("convergence_control" if method=="coulomb" else "deliberately_finite_boundary_control_not_production_accepted"),
    source_file=str(path.relative_to(ROOT)),source_file_sha256=sha(path),source_samples_sha256=data["source_samples_sha256"],
    r=RADII.tolist(),radial_R=enc(R),radial_Rprime=enc(dR),
    scalar22_equatorial_per_q_epsilon=enc(R*S/mu**3),
    propagating=bool(green.propagating),ZI=enc(zi if green.propagating else 0j),up_coefficient=enc(zi),ZH=enc(zh),
    nominal_asymptotic_flux=flux,
    source_free_current_radii=checkr.tolist(),source_free_actual_orbital_current=actualenergy.tolist(),
    actual_current_over_nominal_FI=(float(actualenergy[0]/flux["infinity"]["orbital_energy"]) if green.propagating else None),
    six_term_series_last_term_relative=float(abs(f6-f5)/abs(f6)),
    built_in_last_term_indicator_not_valid_for_forced_order=green.series_last_term_relative,
    wronskian_relative_spread=float(np.max(abs(green.wronskian(np.geomspace(green.rmin,320,100))/green.w0-1))),
    continuous_vs_gauss_zi_relative=float(abs(zi-gausszi)/max(abs(zi),1e-300)),
    continuous_vs_gauss_zh_relative=float(abs(zh-gausszh)/max(abs(zh),1e-300)),
    elapsed_seconds=time.perf_counter()-started)

def main():
 from environment_wake import EnvironmentalWake
 from source_provenance import source_fingerprint,validate_saved_samples
 fingerprint=source_fingerprint()
 paths=[ROOT/"src/environment_radial.py",ROOT/"src/environment_response.py",Path(__file__)]
 inputs={str(p.relative_to(ROOT)):sha(p) for p in paths}
 with ProcessPoolExecutor(4) as pool:
  rows=list(pool.map(calculate,[(orbit,method,order,outer) for orbit in [41.1,42.1] for method,order,outer in CONFIGS]))
 alignmentpath=OUT/"li_field_alignment.json";alignment=json.loads(alignmentpath.read_text())
 fluxpath=OUT/"li_new_radius_flux.json";fluxreference=json.loads(fluxpath.read_text())
 polarpath=ROOT/"docs/environment_reproduction/li_field_polar_reference_20260917.npz"
 polar=np.load(polarpath)
 for p in [alignmentpath,fluxpath,polarpath]:inputs[str(p.relative_to(ROOT))]=sha(p)
 baseline_fields={}
 for orbit in [41.1,42.1]:
  reports=[]
  for relative,h in alignment["inputs_sha256"].items():
   p=ROOT/relative
   if p.suffix!=".json":continue
   data=json.loads(p.read_text());param=data.get("parameters",{})
   if param.get("metric",{}).get("orbital_radius")!=orbit or "samples" not in data:continue
   if sha(p)!=h:raise ValueError("Alignment source artifact changed")
   validate_saved_samples(data,fingerprint);reports.append(data);inputs[relative]=h
  assert len(reports)==18
  wake=EnvironmentalWake(reports,ellmax=5,allow_mixed_discretization=True)
  field=wake.evaluate(RADII[:,None],np.pi/2,polar["phi"][None,:])/(.3**3)
  baseline_fields[orbit]=field
  base=next(r for r in rows if r["orbit"]==orbit and r["boundary_method"]=="coulomb" and r["outer_radius"]==32000)
  reference=next(r for r in fluxreference["rows"] if r["orbit"]==orbit)
  old22=next(m["orbital_energy_infinity"] for m in reference["modes"] if (m["ell"],m["m"])==(2,2))
  nominal_total_baseline=reference["local_infinity_flux"]
  li_flux=reference["references"]["Li"]["reference"]
  refindex=int(np.flatnonzero(polar["rp"]==orbit)[0])
  for row in [r for r in rows if r["orbit"]==orbit]:
   z=dec(row["scalar22_equatorial_per_q_epsilon"]);zb=dec(base["scalar22_equatorial_per_q_epsilon"])
   row["scalar22_amplitude_ratio_to_coulomb32k"]=(abs(z)/abs(zb)).tolist()
   row["scalar22_complex_relative_difference_to_coulomb32k"]=(abs(z-zb)/abs(zb)).tolist()
   changed=field+(z-zb)[:,None]*np.exp(2j*polar["phi"])[None,:]
   row["full18_only_scalar22_replaced"]=[]
   for i,radius in enumerate(RADII):
    indices=np.flatnonzero(np.isclose(polar["radii"],radius,rtol=0,atol=1e-10))
    if not len(indices):
     row["full18_only_scalar22_replaced"].append(dict(radius=float(radius),local_angular_rms=float(np.sqrt(np.mean(abs(changed[i])**2))),reference_available=False))
     continue
    ri=int(indices[0])
    ref=polar["nominal_abs_field"][refindex,ri]
    low=polar["lower_color_interval"][refindex,ri];high=polar["upper_color_interval"][refindex,ri]
    amp=abs(changed[i]);rrms=float(np.sqrt(np.mean(ref**2)))
    row["full18_only_scalar22_replaced"].append(dict(radius=float(radius),
      local_angular_rms=float(np.sqrt(np.mean(amp**2))),reference_nominal_angular_rms=rrms,
      angular_rms_ratio=float(np.sqrt(np.mean(amp**2))/rrms),
      normalized_amplitude_rms_difference=float(np.sqrt(np.mean((amp-ref)**2))/rrms),
      fraction_inside_raster_color_interval=float(np.mean((amp>=low)&(amp<=high)))))
   fi=row["nominal_asymptotic_flux"]["infinity"]["orbital_energy"]
   nominaltotal=nominal_total_baseline-old22+fi
   actualtotal=nominal_total_baseline-old22+row["source_free_actual_orbital_current"][0]
   row["total_infinity_flux_only_scalar22_replaced"]=dict(
     unchanged_other_channels=nominal_total_baseline-old22,
     nominal_asymptotic_total=nominaltotal,actual_current_total=actualtotal,
     Li_vector_plot_reference=li_flux,
     nominal_total_over_Li=nominaltotal/li_flux,actual_current_total_over_Li=actualtotal/li_flux,
     interpretation="For deliberately unconverged boundaries, nominal unit-asymptotic normalization need not equal conserved radial current; neither is a certified retarded infinity result.")
 if source_fingerprint()!=fingerprint:raise ValueError("Production source changed")
 result=dict(status="fixed_fresh_source_outer_boundary_control_completed_not_author_algorithm_identification",
  rows=rows,source_fixed=True,production_files_modified=False,source_provenance=fingerprint,inputs_sha256=inputs,
  source_and_full_field_background="Stationary synchronized a=.8771530275949366; no fixed-a=.88 source replacement.",
  forced_order_note="Only environment_radial.infinity_series is temporarily overridden inside isolated worker processes; restored after construction. This bypasses the production unresolved-series rejection intentionally.",
  leading_boundary="R_up=exp(i k r) r^(-1+i beta/k) times the existing constant phase, with no inverse-r corrections; bound Up may be normalized to one.",
  limitations=["Only scalar22 propagator is changed; all source samples and the other 17 full-field channels are frozen.",
    "Series0/6 finite boundaries are deliberate controls and are not selected or calibrated against the paper.",
    "Closeness of a wrong boundary to reference pixels would not establish the authors' unknown outer radius or algorithm.",
    "Small Wronskian spread establishes integration consistency, not the physical adequacy of the boundary.",
    "Li reference field is quantized original-image color data and flux is a digitized vector curve.",
    "No claim that all retarded-boundary errors are covered, or that this is a full parameter-aligned source calculation."])
 out=OUT/"scalar22_outer_boundary_control.json";out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
 for row in rows:
  print(row["orbit"],row["boundary_method"],row["inverse_r_order"],row["outer_radius"],
    "A ratios",row["scalar22_amplitude_ratio_to_coulomb32k"],
    "nominal/current FI ratio",row["actual_current_over_nominal_FI"],
    "total/Li",row["total_infinity_flux_only_scalar22_replaced"]["nominal_total_over_Li"],
    "ring errors",[v.get("normalized_amplitude_rms_difference") for v in row["full18_only_scalar22_replaced"]],flush=True)
 print(out)
if __name__=="__main__":main()
