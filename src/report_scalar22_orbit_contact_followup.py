"""Scalar22 actual orbit-contact audit at the two near-threshold orbits."""
import concurrent.futures,hashlib,json,time
from pathlib import Path
import numpy as np
from report_orbit_source_matching_contact import work,enc,dec
from environment_source import ThresholdCloud,angular_mode
from environment_radial import RadialGreen
from environment_response import SampledResponse
from lorenz_ghp import KerrGHP
from source_provenance import source_fingerprint,local_dependency_hashes,validate_saved_samples
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'docs/root_cause_followup_20260917'
def main():
 OUT.mkdir(exist_ok=True);cloud=ThresholdCloud(alpha=.3)
 for r0 in [41.1,42.1]:
  base=ROOT/f'docs/field_alignment_20260917/nonstatic/rp{str(r0).replace(".","p")}_sl2_sm2_L6_q12_nr8_h32.json'
  baseline=json.loads(base.read_text());validate_saved_samples(baseline,source_fingerprint())
  config=dict(r0=r0,a=cloud.a,alpha=cloud.mu,quadrature=16,epsilon=.0005,jet_order=8,ellmax=8)
  output=OUT/f'orbit_contact_scalar22_rp{str(r0).replace(".","p")}.json'
  metadata=dict(config=config,baseline=str(base.relative_to(ROOT)),baseline_sha256=hashlib.sha256(base.read_bytes()).hexdigest(),source_provenance=source_fingerprint(),implementation_sha256=local_dependency_hashes(ROOT/'src',[Path(__file__).stem]))
  if output.exists():
   result=json.loads(output.read_text());assert result['metadata']==metadata
  else:result=dict(status='running',metadata=metadata,rows=[])
  def save():
   tmp=output.with_suffix('.tmp');tmp.write_text(json.dumps(result,indent=2)+'\n');tmp.replace(output)
  done={v['ell'] for v in result['rows']};save()
  with concurrent.futures.ProcessPoolExecutor(max_workers=3) as pool:
   futures=[pool.submit(work,(ell,config)) for ell in range(1,9) if ell not in done]
   for future in concurrent.futures.as_completed(futures):
    result['rows'].append(future.result());save();print(r0,'ell',result['rows'][-1]['ell'],'completed',flush=True)
  x,w=np.polynomial.legendre.leggauss(config['quadrature']);tt=np.arccos(x);a=cloud.a
  omega=baseline['parameters']['omega'];S=angular_mode(tt,2,2,a*a*(omega**2-cloud.mu**2))[0]
  R,Rp=cloud.radial(r0)[:,0];Sc,Sp,_=angular_mode(tt,1,1,cloud.c2)
  grad=np.asarray([-1j*cloud.omega*R*Sc,Rp*Sc,R*Sp,1j*R*Sc]).T
  inv=np.asarray([[[z.value for z in row] for row in KerrGHP(r0,t,a,order=2).inv] for t in tt]);sig=r0*r0+a*a*x*x
  g=RadialGreen(a,cloud.mu,omega,2,2,rmax=32000,offset=1e-4,rtol=1e-11,infinity_method='coulomb')
  response=SampledResponse.from_report(baseline);probes=np.array([20.,50.,100.,150.,180.]);original=response.evaluate(probes)[0]
  zh=complex(*baseline['z_h']);zi=complex(*baseline['z_inf']);hz=g.upsol.sol(r0)[0]/g.w0;iz=g.insol.sol(r0)[0]/g.w0
  total=np.zeros((2,len(x),4,4),complex);total1=total.copy();bulk=np.zeros((2,len(x),4),complex);summary=[]
  def project(h):
   jump=h[1]-h[0];trace=np.einsum('kab,kab->k',inv,jump)
   raised=np.einsum('kai,kij,kjb->kab',inv,jump,inv)-.5*inv*trace[:,None,None]
   return 2*np.pi*np.dot(w*S,sig*np.einsum('kb,kb->k',raised[:,1,:],grad))
  for row in sorted(result['rows'],key=lambda x:x['ell']):
   # Shared worker returns conjugate (+mg) data for its historical -mg audit.
   # Undo that conjugation to audit the present +mg1 scalar22 source.
   total+=dec(row['metric_limits']).conjugate();total1+=dec(row['metric_limits_linear']).conjugate();bulk+=dec(row['lorenz_constraint_covariant_at_samples']).conjugate()
   c=project(total);c1=project(total1)
   df=c*np.where(probes<r0,g.insol.sol(probes)[0]*hz,g.upsol.sol(probes)[0]*iz)
   bulkproj=[2*np.pi*np.dot(w*S,sig*np.einsum('kb,kbc,kc->k',v,inv,grad)) for v in bulk]
   summary.append(dict(metric_ellmax=row['ell'],contact_source=enc(c),linear_extrapolated_contact=enc(c1),delta_ZH=enc(c*hz),relative_delta_ZH=float(abs(c*hz/zh)),delta_ZI=enc(c*iz),relative_delta_ZI=float(abs(c*iz/zi)) if zi else None,r=probes.tolist(),relative_complex_field_correction=enc(df/original),relative_field_correction_modulus=abs(df/original).tolist(),bulk_lorenz_source_projection=enc(bulkproj),interpretation='Diagnostic omitted Lorenz contact; not added to physical response. Exact Lorenz matching should make it zero.'))
  result.update(status='completed',summary=summary,limitations=['Finite metric ell cutoff and one Taylor offset; not a global Einstein/source matching proof.','Full baseline uses ell6; ell7/8 rows test contact convergence only, not a refined full scalar response.','Bulk Lorenz samples are local orbit checks, not bounds over all radii.'])
  save();print('SUMMARY',r0,json.dumps(summary[-1]),flush=True)
if __name__=='__main__':main()
