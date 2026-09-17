"""Bounded validation of the new complex Kerr cloud; no environment field run."""
import hashlib,json
from pathlib import Path
import numpy as np
from environment_kerr_quasibound_cloud import GeneralKerrQuasiboundCloud,ComplexScalarAngular,project_complex_source
from environment_source import ThresholdCloud,project_source,kerr_metric
from environment_radial import horizon_series
from environment_schwarzschild_cloud import SchwarzschildCloud
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"docs/root_cause_followup_20260917"
def enc(z):
 z=np.asarray(z,complex);return np.stack([z.real,z.imag],axis=-1).tolist()
def relative(a,b):return float(np.linalg.norm(np.asarray(a)-b)/max(np.linalg.norm(b),1e-300))
def direct_ks_density(cloud,r,theta):
 F,Fp=cloud.ks_radial(r)[:,0];S,Sp=cloud.angular_eigenfunction.evaluate(theta);w=cloud.spectral_omega;a=cloud.a
 phi=F*S;gradient=np.array([-1j*w*phi,Fp*S,F*Sp,1j*cloud.m*phi])
 sigma=r*r+a*a*np.cos(theta)**2;delta=(r-cloud.rp)*(r-cloud.rm)
 inverse=np.array([[-1-2*r/sigma,2*r/sigma,0,0],[2*r/sigma,delta/sigma,0,a/sigma],[0,0,1/sigma,0],[0,a/sigma,0,1/(sigma*np.sin(theta)**2)]])
 metric=np.linalg.inv(inverse)
 L=np.conjugate(gradient)@inverse@gradient+cloud.mu**2*abs(phi)**2
 stress=np.outer(gradient.conjugate(),gradient)+np.outer(gradient,gradient.conjugate())-metric*L
 energy=-sigma*(inverse@stress)[0,0].real
 raised=inverse@gradient;current=-1j*(phi.conjugate()*raised-phi*raised.conjugate())
 charge=sigma*current[0].real
 fn,sn=abs(F)**2,abs(S)**2;cross=np.imag(F.conjugate()*Fp)
 analyticE=((sigma+2*r)*abs(w)**2+cloud.mu**2*sigma+cloud.m**2/np.sin(theta)**2)*fn*sn+delta*abs(Fp)**2*sn+fn*abs(Sp)**2+2*a*cloud.m*cross*sn
 analyticQ=2*((sigma+2*r)*w.real*fn+2*r*cross)*sn
 return dict(r=r,theta=theta,direct_stress_energy=energy,integrand_energy=float(analyticE),direct_Noether_charge=charge,integrand_charge=float(analyticQ),energy_relative=abs(energy-analyticE)/abs(energy),charge_relative=abs(charge-analyticQ)/abs(charge))
def checks(cloud):
 samples=[]
 for r in [cloud.rp+.0005,3.,20.,100.]:
  theta=1.1;field,H,inv=cloud.hessian(r,theta)
  kg=np.einsum("ij,ij->",inv,H)-cloud.mu**2*field
  geom,jet=cloud.jet(r,theta,order=4,return_geometry=True)
  jetkg=geom.scalar_wave(jet).value-cloud.mu**2*jet.value
  spectral=cloud.spectral_omega;temporal=complex(cloud.omega)
  expected=(inv[0,0]*(spectral**2-temporal**2)+2*inv[0,3]*cloud.m*(temporal-spectral))*field
  R,Rp,Rpp=cloud.radial_state(r)[:,0]
  step=min(.001*(r-cloud.rp),.005*max(r,1))
  offsets=np.array([-2,-1,1,2])*step+r
  diff=np.dot([1,-8,8,-1],cloud.radial(offsets)[1])/(12*step)
  delta,dp,V=__import__("environment_cloud").radial_coefficients(r,cloud.a,cloud.mu,spectral,1,cloud.lam)
  residual=delta*diff+dp*Rp+V*R
  scale=abs(delta*diff)+abs(dp*Rp)+abs(V*R)
  samples.append(dict(r=r,field=enc(field),kg_residual=enc(kg),kg_relative_to_mu2_field=float(abs(kg)/(cloud.mu**2*abs(field))),
   expected_freezing_defect=enc(expected),kg_minus_predicted_defect_relative=float(abs(kg-expected)/(cloud.mu**2*abs(field))),
   jet_kg_relative=float(abs(jetkg)/(cloud.mu**2*abs(field))),
   jet_minus_hessian_kg_relative=float(abs(jetkg-kg)/(cloud.mu**2*abs(field))),
   jet_field_relative=relative(jet.value,field),jet_Rtheta_derivative_relative=relative(jet.derivative(0).value,Rp*cloud.angular_eigenfunction.evaluate(theta)[0]),
   independent_five_point_radial_residual_relative=float(abs(residual)/scale)))
 return samples
def main():
 OUT.mkdir(exist_ok=True)
 exact=GeneralKerrQuasiboundCloud(rtol=2e-13,outer_efolds=75.,radial_points=6001)
 coarse=GeneralKerrQuasiboundCloud(rtol=2e-11,offset=1e-5,outer_efolds=60.,radial_points=3001,horizon_order=4,infinity_order=6)
 frozen=GeneralKerrQuasiboundCloud(freeze_growth=True,rtol=2e-13,outer_efolds=75.,radial_points=6001)
 rr=[exact.rp+.001,3.,20.,100.,320.]
 result=dict(status="bounded_general_Kerr_quasibound_cloud_validation",exact_cloud_provenance=exact.provenance,
  matching_residual=enc(exact.matching_log_derivative_residual),
  independent_mass_reintegration=dict(grid6001=exact.ks_integrals(6001),grid12001=exact.ks_integrals(12001)),
  horizon_balance=exact.horizon_balance(12001),
  endpoint_tolerance_refinement=dict(coarse_provenance=coarse.provenance,samples=[dict(r=r,radial_state_relative=relative(coarse.radial_state(r),exact.radial_state(r))) for r in rr]),
  direct_KS_stress_tensor_checks=[direct_ks_density(exact,r,t) for r,t in [(3.,.7),(20.,1.1),(100.,2.)]],
  exact_local_checks=checks(exact),frozen_local_checks=checks(frozen),
  freeze_mass_defect=dict(spectral_energy=frozen.ks_integrals(12001)[0],temporally_frozen_energy_on_same_slice=frozen.ks_integrals(12001,frequency=frozen.omega)[0]),
  limitations=["This validates the cloud and local source interfaces only; no new metric/source grid or full forced field was calculated.",
    "Exact complex-frequency source projection uses the bilinear angular dual; real-frequency RadialGreen/mode_flux remain inapplicable without explicit approximation.",
    "freeze_growth is temporal-only AFTER the complex spatial profile and spectral KS mass normalization; it has a measured KG defect and is not a synchronized or exact stationary solution."])
 R=exact.radial(exact.rmin)[:,0]/exact.amplitude
 original=np.array(horizon_series(exact.rmin,exact.a,exact.mu,exact.spectral_omega,1,exact.lam,order=6))
 result["independent_BL_horizon_series_comparison"]=dict(radius=exact.rmin,relative_state_difference=relative(R,original),new=enc(R),existing_independent_series=enc(original))
 old=ThresholdCloud(alpha=.3,radial_points=6001)
 sync=GeneralKerrQuasiboundCloud(a=old.a,omega=old.omega,radial_points=6001,rtol=2e-13)
 result["synchronous_limit"]=dict(a=old.a,omega=old.omega,KS_energy=sync.ks_integrals(12001)[0],old_BL_energy=old.integrals(12001)[0],
  samples=[dict(r=r,radial_state_relative=relative(sync.radial(r),old.radial(r)),hessian_relative=relative(sync.hessian(r,1.1)[1],old.hessian(r,1.1)[1])) for r in [old.rp+.001,3.,20.,100.,320.]])
 oldsch=SchwarzschildCloud(alpha=.3,freeze_decay=False,radial_points=6001)
 sch=GeneralKerrQuasiboundCloud(a=0.,omega=oldsch.spectral_omega,radial_points=6001,rtol=2e-13)
 result["Schwarzschild_limit"]=dict(omega=enc(oldsch.spectral_omega),new_balance=sch.horizon_balance(12001),old_balance=oldsch.decay_balance(),
  samples=[dict(r=r,BL_radial_relative=relative(sch.radial(r),oldsch.radial(r)),KS_radial_relative=relative(sch.ks_radial(r),np.array(oldsch.ks_radial(r)))) for r in [2.001,3.,20.,100.,320.]])
 x,weights=np.polynomial.legendre.leggauss(80);theta=np.arccos(x)
 A1=ComplexScalarAngular(1,1,.1+.2j);A3=ComplexScalarAngular(3,1,.1+.2j)
 s1=A1.evaluate(theta)[0];s3=A3.evaluate(theta)[0]
 result["complex_angular_dual_test"]=dict(test_c2=enc(.1+.2j),physical_cloud_c2=enc(exact.c2),
  hermitian_norm=float((2*np.pi*np.dot(weights,abs(s1)**2)).real),
  dual_self=enc(2*np.pi*np.dot(weights,A1.dual(theta)*s1)),
  dual_other=enc(2*np.pi*np.dot(weights,A1.dual(theta)*s3)),
  incorrect_Hermitian_cross_projection=enc(2*np.pi*np.dot(weights,s1.conjugate()*s3)))
 def toy(r,t):return np.diag([.01/r,.02/r,.03*r,.04*r*np.sin(t)**2]).astype(complex)
 w1,J1=project_source(frozen,[3.,20.],42.1,2,2,toy,ntheta=24)
 w2,J2=project_complex_source(frozen,[3.,20.],42.1,2,2,toy,ntheta=24)
 we,Je=project_complex_source(exact,[3.,20.],42.1,2,2,toy,ntheta=24)
 result["source_interface_test"]=dict(metric="Fixed arbitrary symmetric toy tensor only; no Lorenz/gauge or physical source claim.",
  frozen_existing_projection=enc(J1),frozen_explicit_projection=enc(J2),relative_difference=relative(J1,J2),
  exact_complex_omega=enc(we),exact_complex_projection=enc(Je))
 result["implementation_sha256"]={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/"src/environment_kerr_quasibound_cloud.py",ROOT/"src/environment_source.py",ROOT/"src/environment_cloud.py",ROOT/"src/environment_radial.py",ROOT/"src/environment_schwarzschild_cloud.py",ROOT/"src/lorenz_mode_jet.py"]}
 out=OUT/"general_kerr_quasibound_cloud_validation.json";out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
 print(json.dumps({k:v for k,v in result.items() if k not in ["implementation_sha256","exact_cloud_provenance","source_interface_test"]},indent=2))
if __name__=="__main__":main()
