"""Only ell_g=1 full-current spin1 ablation; never modify the physical solver.

Full and vector sources share saved radial nodes, weights, cloud and massive
Green kernel. Negative metric m is obtained by conjugating only the positive-m
BL metric. The cloud is NOT conjugated and there is no factor-of-two adjustment.
A deleted gauge sector need not remain Lorenz or solve the sourced Einstein
equation; 'legacy' means a diagnostic of a proposed historical implementation.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from environment_source import ThresholdCloud,project_source
from environment_radial import RadialGreen
from environment_cloud import mode_flux
from lorenz_metric import spin1_metric
from source_provenance import source_fingerprint,validate_saved_samples,local_dependency_hashes
ROOT=Path(__file__).resolve().parents[1]
def enc(z):return [float(np.real(z)),float(np.imag(z))]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
class VectorDipole:
    def __init__(self,r0,a,order):self.r0,self.a,self.order=r0,a,order
    def __call__(self,r,t):
        _,h=spin1_metric(r,t,self.r0,self.a,1,1,order=self.order,full_current=True)
        return np.array([[v.value for v in row] for row in h],complex).conjugate()

def main():
    p=argparse.ArgumentParser();p.add_argument('--baseline',default='docs/environment_reproduction/fresh_20260915_L18_mg-1_sl0.json');p.add_argument('--output',required=True)
    p.add_argument('--angular-order',type=int);p.add_argument('--jet-order',type=int,default=6);p.add_argument('--green-rtol',type=float,default=1e-11)
    args=p.parse_args();baseline=Path(args.baseline);old=json.loads(baseline.read_text());par=old['parameters'];orbit=par['metric']['orbital_radius'];m=par['metric']['m_g']
    if (par['alpha'],orbit,m,par['scalar_ell'],par['scalar_m'])!=(.3,20.,-1,0,0):raise ValueError('This first experiment is alpha.3/rp20/mg-1/scalar00')
    validate_saved_samples(old,source_fingerprint())
    out=Path(args.output)
    if out.exists():raise FileExistsError('Preserve existing ablation result; use a new output name')
    cloud=ThresholdCloud(alpha=par['alpha']);nq=args.angular_order or par['angular_order'];radii=np.array([x['r'] for x in old['samples']]);weights=np.array([x['weight'] for x in old['samples']]);Jfull=np.array([complex(*x['source']) for x in old['samples']])
    if cloud.a!=par['metric']['a'] or cloud.mass!=par['cloud_mass']:raise ValueError('Baseline cloud normalization/spin mismatch')
    omega=cloud.omega-1/(orbit**1.5+cloud.a)
    if omega!=par['omega']:raise ValueError('Frequency differs from baseline')
    green=RadialGreen(cloud.a,cloud.mu,omega,0,0,rmax=par['green_outer_radius'],offset=par['green_horizon_offset'],rtol=args.green_rtol,infinity_method=par.get('infinity_method','series'))
    kernel=green.upsol.sol(radii)[0]/green.w0;Zfull=np.dot(weights,kernel*Jfull)
    relative=abs(Zfull-complex(*old['z_h']))/abs(Zfull)
    if relative>2e-9:raise ValueError('Fresh Green integration fails stored full-source amplitude check')
    metric=VectorDipole(orbit,cloud.a,args.jet_order);start=time.perf_counter()
    result=dict(status='sampling_vector_only',parameters=par,controls=dict(vector_angular_order=nq,vector_jet_order=args.jet_order,green_rtol=args.green_rtol),baseline=str(baseline),baseline_sha256=sha(baseline),source_provenance=source_fingerprint(),implementation_sha256=local_dependency_hashes(ROOT/'src',['report_legacy_vector_ablation']),baseline_Z_reintegration_relative=relative,operation='h_legacy=h_physical-conjugate(h_spin1(ell1,m+1,full_current=True)); cloud,trace,kappa,higherell,Green,prefactor unchanged',samples=[],limitations=['This tests a proposed legacy missing vector implementation; deleting a gauge sector need not retain Lorenz gauge or sourced Einstein matching.','Same physical full-source radial grid and weights. If vector angular order differs, only the vector subtraction is angular-refined; full baseline source remains its recorded angular order.','Scalar00 channel only, not the total horizon or infinity sum. No claim of provenance of the2025 figure follows from an ablation result.'])
    def save():
        tmp=out.with_suffix('.tmp');tmp.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');tmp.replace(out)
    save()
    for i,r in enumerate(radii):
        frequency,values=project_source(cloud,[r],orbit,0,0,metric,ntheta=nq)
        if frequency!=omega or not np.isfinite(values).all():raise ValueError('Bad vector source')
        result['samples'].append(dict(r=float(r),weight=float(weights[i]),source_full=enc(Jfull[i]),source_vector=enc(values[0]),source_legacy=enc(Jfull[i]-values[0]),green_up_over_W=enc(kernel[i])))
        if (i+1)%8==0 or i+1==len(radii):save();print(f'vector {i+1}/{len(radii)} r={r:.8g}',flush=True)
    JV=np.array([complex(*s['source_vector']) for s in result['samples']]);ZV=np.dot(weights,kernel*JV);Zlegacy=np.dot(weights,kernel*(Jfull-JV))
    flux=lambda z:mode_flux(omega,0,cloud.omega,cloud.m,cloud.mu,cloud.a,0j,z)['horizon']['orbital_energy']
    constant=flux(1+0j);FF,FV,FL=flux(Zfull),flux(ZV),flux(Zlegacy);interference=2*constant*np.real(Zlegacy*np.conjugate(ZV))
    result.update(status='completed_same_grid_vector_ablation',z_h_full=enc(Zfull),z_h_vector=enc(ZV),z_h_legacy=enc(Zlegacy),flux_h_full=FF,flux_h_vector=FV,flux_h_legacy=FL,interference_full_legacy_vector=interference,flux_identity_residual=FF-FL-FV-interference,amplitude_identity_residual=enc(Zfull-ZV-Zlegacy),legacy_over_physical_flux=FL/FF,relative_flux_change=(FL-FF)/abs(FF),vector_amplitude_over_full=enc(ZV/Zfull),flux_prefactor=constant,source_integral_cancellation=dict(full=float(np.sum(abs(weights*kernel*Jfull))/abs(Zfull)),vector=float(np.sum(abs(weights*kernel*JV))/abs(ZV)),legacy=float(np.sum(abs(weights*kernel*(Jfull-JV)))/abs(Zlegacy))),wronskian_relative_spread=float(np.max(abs(green.wronskian(radii)/green.w0-1))),elapsed_seconds=time.perf_counter()-start)
    save();print(json.dumps({k:result[k] for k in ['z_h_full','z_h_vector','z_h_legacy','flux_h_full','flux_h_legacy','legacy_over_physical_flux','interference_full_legacy_vector']},indent=2),flush=True)
if __name__=='__main__':main()
