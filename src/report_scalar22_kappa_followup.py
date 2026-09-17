"""Actual scalar22 source correction from independent analytic kappa variation."""
from pathlib import Path
import hashlib,json,time
import numpy as np
from environment_source import ThresholdCloud,angular_mode,connection
from environment_response import SampledResponse
from environment_trace_variation import TraceMassVariation
from lorenz_kappa import kappa_jet
from lorenz_ghp import KerrGHP
from report_dipole_kappa_flux_audit import scalar_hessian
from source_provenance import source_fingerprint,local_dependency_hashes,validate_saved_samples
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'docs/root_cause_followup_20260917'
def enc(z):
 z=np.asarray(z,complex);return np.stack([z.real,z.imag],axis=-1).tolist()
def main():
 OUT.mkdir(exist_ok=True);cloud=ThresholdCloud(alpha=.3);x,w=np.polynomial.legendre.leggauss(12);theta=np.arccos(x)
 cases=[(1,2000.),(1,8000.),(3,2000.),(3,8000.)]
 for rp in [41.1,42.1]:
  baseline=ROOT/f'docs/field_alignment_20260917/nonstatic/rp{str(rp).replace(".","p")}_sl2_sm2_L6_q12_nr8_h32.json'
  d=json.loads(baseline.read_text());validate_saved_samples(d,source_fingerprint());p=d['parameters'];wg=1/(rp**1.5+cloud.a)
  s=angular_mode(theta,2,2,cloud.a**2*(p['omega']**2-cloud.mu**2))[0]
  r=np.asarray([v['r'] for v in d['samples']]);J=np.asarray([complex(*v['source']) for v in d['samples']]);weights=np.asarray([v['weight'] for v in d['samples']])
  variations=[TraceMassVariation(rp,cloud.a,ell,1,rmax=outer) for ell,outer in cases]
  correction=np.zeros((len(cases),len(r)),complex);st=time.monotonic()
  output=OUT/f'kappa_scalar22_rp{str(rp).replace(".","p")}.json'
  for i,rr in enumerate(r):
   for t,xx,ww,ss in zip(theta,x,w,s):
    g=KerrGHP(float(rr),float(t),cloud.a,omega=wg,m=1,order=2)
    old={ell:kappa_jet(g,rp,ell) for ell in (1,3)}
    inv,gamma=connection(rr,t,cloud.a);hc=cloud.hessian(rr,t)[1];raised=inv@hc@inv
    factor=2*np.pi*ww*ss*(rr*rr+cloud.a**2*xx*xx)
    for j,((ell,outer),v) in enumerate(zip(cases,variations)):
     dh=-2j/wg*scalar_hessian(v.kappa_jet(g)-old[ell],wg,1,gamma)
     correction[j,i]+=factor*np.einsum('ab,ab->',dh,raised)
   if (i+1)%16==0:print(rp,i+1,'/',len(r),flush=True)
  original=SampledResponse.from_report(d);green=original.green;probes=np.array([20.,50.,70.,100.,150.,180.]);R0=original.evaluate(probes)[0]
  rows=[]
  for j,(ell,outer) in enumerate(cases):
   dn=SampledResponse(green,p['source_panels'],r,correction[j],log_first=True);delta=dn.evaluate(probes)[0]
   z=original.up_coefficient;zh=original.horizon_coefficient
   rows.append(dict(ell=ell,analytic_trace_outer_radius=outer,delta_source=enc(correction[j]),delta_ZH=enc(dn.horizon_coefficient),delta_ZI=enc(dn.up_coefficient),relative_delta_ZH=float(abs(dn.horizon_coefficient/zh)),relative_delta_up_coefficient=float(abs(dn.up_coefficient/z)),infinity_flux_relative_change=float(abs(1+dn.up_coefficient/z)**2-1) if green.propagating else None,r=probes.tolist(),relative_complex_field_correction=enc(delta/R0),relative_field_correction_modulus=abs(delta/R0).tolist()))
  result=dict(status='completed_isolated_spin0_kappa_source_control',baseline=str(baseline.relative_to(ROOT)),baseline_sha256=hashlib.sha256(baseline.read_bytes()).hexdigest(),source_provenance=source_fingerprint(),implementation_sha256=local_dependency_hashes(ROOT/'src',[Path(__file__).stem]),rows=rows,elapsed_seconds=time.monotonic()-st,limitations=['Only kappa at metric ell1 and ell3 is varied; their corrections have not been inserted into production.','Original source quadrature/angular cutoff retained.','Source correction includes tensor Hessians explicitly; no fitted field amplitude.'])
  output.write_text(json.dumps(result,indent=2)+'\n');print('DONE',rp,[(v['ell'],v['analytic_trace_outer_radius'],v['infinity_flux_relative_change'],max(v['relative_field_correction_modulus'])) for v in rows],flush=True)
if __name__=='__main__':main()
