"""Independent scalar angular/phase audit for Fig.7, without refitting the field.

Read-only with respect to production modules and saved mode responses.
"""
import argparse
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy import special
from environment_source import angular_mode, project_source
from source_provenance import source_fingerprint

ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'docs/environment_reproduction'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def pair(z):return [float(np.real(z)),float(np.imag(z))]
def scaled_error(a,b):return float(np.max(np.abs(a-b))/max(float(np.max(np.abs(b))),1e-300))
def standard(theta,ell,m):return special.sph_harm(m,ell,0.,theta).real

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=FOLDER/'particle_phase_audit_20260916.json');args=parser.parse_args()
 theta=np.r_[np.linspace(.035,np.pi-.035,31),np.pi/2];x,w=np.polynomial.legendre.leggauss(80);qt=np.arccos(x)
 spherical=[]
 for ell in range(19):
  for m in range(-ell,ell+1):
   s,ds,lam=angular_mode(theta,ell,m,0.);ref=standard(theta,ell,m)
   dr=np.zeros_like(theta)
   if m<ell:dr+=.5*np.sqrt((ell-m)*(ell+m+1))*standard(theta,ell,m+1)
   if m>-ell:dr-=.5*np.sqrt((ell+m)*(ell-m+1))*standard(theta,ell,m-1)
   sq=angular_mode(qt,ell,m,0.)[0]
   pos=angular_mode(theta,ell,abs(m),0.)[0]
   expected=pos*((-1)**abs(m) if m<0 else 1.)
   spherical.append(dict(ell=ell,m=m,value_relative_max=scaled_error(s,ref),derivative_relative_max=scaled_error(ds,dr) if ell else float(np.max(abs(ds))),normalization_error=abs(float(2*np.pi*np.dot(w,sq*sq))-1.),negative_m_relation_error=float(np.max(abs(s-expected))),eigenvalue_error=float(abs(lam-ell*(ell+1)))))
 spheroidal=[]
 for c2 in (-.07,-.003,.02,.2):
  c=np.sqrt(abs(c2));fun=special.pro_ang1 if c2<0 else special.obl_ang1;cvfun=special.pro_cv if c2<0 else special.obl_cv
  for ell in range(13):
   for m in range(-ell,ell+1):
    # SciPy angular functions are unnormalized; independently fix their
    # spherical-limit sign by a positive overlap with standard Y_lm.
    sq=fun(abs(m),ell,c,x)[0]
    overlap=float(np.dot(w,sq*standard(qt,ell,m)))
    factor=np.sign(overlap)/np.sqrt(2*np.pi*np.dot(w,sq*sq))
    sr,dx=fun(abs(m),ell,c,np.cos(theta));sr=sr*factor;dsr=-np.sin(theta)*dx*factor
    s,ds,lam=angular_mode(theta,ell,m,c2)
    spheroidal.append(dict(ell=ell,m=m,c2=c2,value_relative_max=scaled_error(s,sr),derivative_relative_max=scaled_error(ds,dsr),eigenvalue_absolute_error=float(abs(lam-cvfun(abs(m),ell,c)))))
 # Exercise the actual production projector with an independent manufactured
 # spherical source. This tests negative-m signs and 2pi normalization.
 fake=SimpleNamespace(a=0.,omega=.29629324847975716,mu=.3,m=1,lorenz_source=lambda r,t,h:h)
 projections=[]
 for ell in range(19):
  for m in range(-ell,ell+1):
   coeff=(1+.07*m)+1j*(.2+.03*ell)
   other=ell+1
   def metric(r,t):return (coeff*standard(t,ell,m)+(.23-.31j)*standard(t,other,m))/(r*r)
   _,out=project_source(fake,[7.],20.,ell,m,metric,ntheta=80)
   projections.append(dict(ell=ell,m=m,expected=pair(coeff),actual=pair(out[0]),relative_error=float(abs(out[0]-coeff)/abs(coeff))))
 particle_path=FOLDER/'particle_field_L18.json';markers_path=FOLDER/'paper_figure7_markers.json';spherical_path=FOLDER/'particle_spherical_projection.json'
 particle=json.loads(particle_path.read_text());markers={v['ell']:v['plotted_value'] for v in json.loads(markers_path.read_text())['markers']};sphere=json.loads(spherical_path.read_text());sphere_rows={v['ell']:v for v in sphere['rows']}
 physical=[];input_checks=[]
 for shell in particle['multipoles']:
  if not shell['complete']:continue
  ell=shell['ell'];radial=0j;field=0j;signchanged_radial=0j;rephased_field=0j;angular_errors=[]
  for comp in shell['components']:
   path=FOLDER/comp['file'];data=json.loads(path.read_text());par=data['parameters'];m=comp['m'];a=par['metric']['a'];alpha=par['alpha'];omega=par['omega']
   z=complex(*comp['radial']);s=float(angular_mode(np.pi/2,ell,m,a*a*(omega*omega-alpha*alpha))[0]);phase=(-1)**((ell-m)//2)
   input_checks.append(dict(file=comp['file'],sha256=sha(path),matches_particle_hash=sha(path)==particle['inputs_sha256'][comp['file']],source_provenance=data.get('source_provenance',{}).get('sha256')))
   angular_errors.append(abs(s-comp['equatorial_angular_value']));radial+=z;field+=z*s;signchanged_radial+=phase*z;rephased_field+=(phase*z)*(phase*s)
  scale=alpha**3;mark=markers[ell];rec=dict(ell=ell,radial_sum=pair(radial),physical_field=pair(field),physical_per_epsilon_q=float(abs(field)/scale),radial_per_epsilon_q=float(abs(radial)/scale),paper_marker=mark,physical_over_marker=float(abs(field)/scale/mark),radial_over_marker=float(abs(radial)/scale/mark),max_saved_angular_value_error=max(angular_errors),saved_field_relative_error=float(abs(field-complex(*shell['particle_field_sum']))/abs(field)),consistent_basis_sign_change_field_relative_error=float(abs(rephased_field-field)/abs(field)),same_sign_change_radial_only_sum_ratio=float(abs(signchanged_radial)/abs(radial)))
  if ell in sphere_rows:
   sfield=complex(*sphere_rows[ell]['spherical_field']);rec.update(spherical_field=pair(sfield),spherical_vs_spheroidal_shell_relative=float(abs(sfield-field)/abs(field)))
  physical.append(rec)
 summary=dict(spherical_modes=len(spherical),max_spherical_value_relative=max(v['value_relative_max'] for v in spherical),max_spherical_derivative_relative=max(v['derivative_relative_max'] for v in spherical),max_normalization_error=max(v['normalization_error'] for v in spherical),max_negative_m_relation_error=max(v['negative_m_relation_error'] for v in spherical),spheroidal_modes=len(spheroidal),max_spheroidal_value_relative=max(v['value_relative_max'] for v in spheroidal),max_spheroidal_derivative_relative=max(v['derivative_relative_max'] for v in spheroidal),max_spheroidal_eigenvalue_absolute=max(v['eigenvalue_absolute_error'] for v in spheroidal),projector_modes=len(projections),max_projector_relative_error=max(v['relative_error'] for v in projections),all_input_hashes_match=all(v['matches_particle_hash'] for v in input_checks))
 summary['all_inputs_have_source_provenance']=all(v['source_provenance'] is not None for v in input_checks)
 summary['declared_tolerances_passed']=summary['max_spherical_value_relative']<1e-10 and summary['max_spherical_derivative_relative']<1e-10 and summary['max_normalization_error']<1e-10 and summary['max_spheroidal_value_relative']<1e-9 and summary['max_spheroidal_derivative_relative']<1e-9 and summary['max_projector_relative_error']<1e-10 and summary['all_input_hashes_match']
 pdf=ROOT/'outputs/environment_reference/Flux_l_mode_convergence_particle.pdf'
 summary['marker_source_pdf_hash_matches']=sha(pdf)==json.loads(markers_path.read_text())['sha256']
 result=dict(status='current_angular_phase_tests_pass_historical_particle_definition_comparison_unresolved',current_source_provenance=source_fingerprint(),summary=summary,spherical=spherical,spheroidal=spheroidal,manufactured_source_projection=projections,particle_shells=physical,input_files=input_checks,input_sha256={str(p.relative_to(ROOT)):sha(p) for p in [particle_path,markers_path,spherical_path,ROOT/'src/environment_source.py',ROOT/'src/environment_lorenz_mode.py',ROOT/'src/report_particle_field.py',ROOT/'outputs/environment_reference/main_PRL.tex',pdf]},implementation_sha256=sha(Path(__file__)),conventions=dict(spherical='S_lm(theta,0)=Y_lm(theta,0), including Condon-Shortley phase',normalization='2pi integral_-1^1 |S|^2 dx=1',negative_m='S_l,-m=(-1)^m S_l,m at identical real c2; this is NOT a conjugacy identity for the driven complex-cloud radial fields',projection='J_lm=2pi integral S_lm Sigma source_m dx; S is real in this implementation',particle='sum_m R_lm(rp) S_lm(pi/2) at t=phi=0',literal_paper_definition='main_PRL.tex:748 writes phi_l=sum_m phi_lm, where lines468/718 define phi_lm as radial coefficient'),limitations=['No recomputation of metric, cloud, or radial response in this audit. Saved response provenance is recorded, not upgraded. Every particle response input currently lacks source_provenance; matching file hashes establish unchanged historical inputs, not current-production provenance.','Independent angular tests do not validate the metric reconstruction or finite radial source quadrature.','Literal radial-only sum depends on the chosen angular basis phases; physical reconstructed field does not. Original environment source/plotting code is unavailable.','Figure7 markers are from historical arXiv v1 and have no formal digitization uncertainty.','Spherical comparison uses the already completed finite even-shell reprojection; omitted ell>12 remains omitted.'])
 args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(summary,indent=2));print(json.dumps(physical,indent=2))
if __name__=='__main__':main()
