"""Direct one-sided production metric audit of the orbit delta source.

No physical module is modified. Negative metric m is conjugated, cloud is not.
The delta coefficient is [bar h^{rb}] d_b Phi0, projected with Sigma S_00.
"""
import argparse,concurrent.futures,hashlib,json,os,time
from pathlib import Path
import numpy as np
from environment_source import ThresholdCloud,angular_mode
from environment_radial import RadialGreen
from lorenz_metric import nonstatic_metric
from lorenz_ghp import KerrGHP
from lorenz_tensor import lorenz_constraint
from source_provenance import local_dependency_hashes,source_fingerprint
ROOT=Path(__file__).resolve().parents[1]

def enc(z):
 a=np.asarray(z,complex);return np.stack([a.real,a.imag],axis=-1).tolist()
def dec(z):
 a=np.asarray(z);return a[...,0]+1j*a[...,1]
def work(task):
 ell,c=task;st=time.monotonic();r0=c['r0'];a=c['a'];eps=c['epsilon'];q=c['quadrature'];order=c['jet_order']
 x,weights=np.polynomial.legendre.leggauss(q);theta=np.arccos(x);limits=[];linear=[];bulk=[]
 for sign in (-1,1):
  r=r0+sign*eps;dr=-sign*eps;side=[];sidelinear=[];sidebulk=[]
  for t in theta:
   g,h=nonstatic_metric(r,t,r0,a,ell,1,order=order)
   v=np.asarray([[z.value for z in row] for row in h]);v1=np.asarray([[z.derivative(0).value for z in row] for row in h]);v2=np.asarray([[z.derivative(0).derivative(0).value for z in row] for row in h])
   side.append((v+dr*v1+dr*dr*v2/2).conjugate());sidelinear.append((v+dr*v1).conjugate())
   con=lorenz_constraint(g,h)
   sidebulk.append(np.asarray([z.value for z in con]).conjugate())
  limits.append(side);linear.append(sidelinear);bulk.append(sidebulk)
 return dict(ell=ell,metric_limits=enc(limits),metric_limits_linear=enc(linear),lorenz_constraint_covariant_at_samples=enc(bulk),elapsed_seconds=time.monotonic()-st)

def summarize(rows,config,cloud,baseline):
 x,w=np.polynomial.legendre.leggauss(config['quadrature']);tt=np.arccos(x);r0=config['r0'];a=cloud.a
 omega=cloud.omega-1/(r0**1.5+a);S=angular_mode(tt,0,0,a*a*(omega*omega-cloud.mu*cloud.mu))[0]
 R,Rp=cloud.radial(r0)[:,0];Sb,Sbp,_=angular_mode(tt,1,1,cloud.c2)
 grad=np.array([-1j*cloud.omega*R*Sb,Rp*Sb,R*Sbp,1j*R*Sb]).T
 inv=np.asarray([[[z.value for z in row] for row in KerrGHP(r0,t,a,order=2).inv] for t in tt]);sigma=r0*r0+a*a*x*x
 pars=baseline['parameters'];green=RadialGreen(a,cloud.mu,omega,0,0,rmax=pars['green_outer_radius'],offset=pars['green_horizon_offset'],rtol=1e-11,infinity_method=pars.get('infinity_method','series'))
 kernel=green.upsol.sol(r0)[0]/green.w0;Z=complex(*baseline['z_h']);total=np.zeros((2,len(x),4,4),complex);total_lin=total.copy();out=[]
 def coefficient(h):
  jump=h[1]-h[0];trace=np.einsum('kab,kab->k',inv,jump);raised=np.einsum('kai,kij,kjb->kab',inv,jump,inv)-.5*inv*trace[:,None,None]
  integrand=sigma*np.einsum('kb,kb->k',raised[:,1,:],grad)
  return 2*np.pi*np.dot(w*S,integrand),np.linalg.norm(jump,axis=(1,2)),np.linalg.norm(h,axis=(2,3)).max(axis=0)
 for row in sorted(rows,key=lambda r:r['ell']):
  total+=dec(row['metric_limits']);total_lin+=dec(row['metric_limits_linear']);contact,jumps,scale=coefficient(total);linear,_,_=coefficient(total_lin);dz=kernel*contact
  out.append(dict(ellmax=row['ell'],contact_J=enc(contact),contact_J_linear_extrapolation=enc(linear),delta_z_h_contact=enc(dz),relative_delta_z_h=float(abs(dz/Z)),signed_flux_relative_change=float(abs(1+dz/Z)**2-1),max_metric_jump_norm=float(max(jumps)),max_metric_jump_over_side_norm=float(max(jumps/np.maximum(scale,1e-300))),scalar_weighted_jump_bound=float(2*np.pi*np.dot(w*abs(S),jumps)),taylor_linear_quadratic_contact_difference=float(abs(contact-linear))))
 return dict(omega=omega,green_Up_over_W_at_orbit=enc(kernel),baseline_z_h=enc(Z),cumulative_ell_results=out)

def main():
 p=argparse.ArgumentParser();p.add_argument('--ellmax',type=int,default=18);p.add_argument('--quadrature',type=int,default=18);p.add_argument('--epsilon',type=float,default=1e-3);p.add_argument('--jet-order',type=int,default=8);p.add_argument('--workers',type=int,default=3);p.add_argument('--output',required=True);p.add_argument('--baseline',default='docs/environment_reproduction/fresh_20260915_L18_mg-1_sl0.json');args=p.parse_args()
 out=Path(args.output);base=Path(args.baseline);baseline=json.loads(base.read_text());p0=baseline['parameters'];cloud=ThresholdCloud(alpha=p0['alpha'])
 if (p0['metric']['m_g'],p0['scalar_ell'],p0['scalar_m'],p0['cloud_mass'])!=(-1,0,0,1.):raise ValueError('Only unit-mass mg-1 scalar00 is supported')
 if p0['metric']['a']!=cloud.a:raise ValueError('cloud spin mismatch')
 config=dict(r0=p0['metric']['orbital_radius'],a=cloud.a,alpha=cloud.mu,quadrature=args.quadrature,epsilon=args.epsilon,jet_order=args.jet_order,ellmax=args.ellmax)
 hashes=local_dependency_hashes(ROOT/'src',['report_orbit_source_matching_contact']);metadata=dict(config=config,baseline_sha256=hashlib.sha256(base.read_bytes()).hexdigest(),implementation_sha256=hashes,source_provenance=source_fingerprint())
 if out.exists():
  d=json.loads(out.read_text())
  if d['metadata']!=metadata:raise ValueError('Cannot resume changed implementation/configuration')
 else:d=dict(status='running_direct_one_sided_contact_audit',metadata=metadata,rows=[])
 done={r['ell'] for r in d['rows']};tasks=[(l,config) for l in range(1,args.ellmax+1) if l not in done]
 def save():
  temp=out.with_suffix('.tmp');temp.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n');temp.replace(out)
 save()
 with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
  futures={pool.submit(work,t):t[0] for t in tasks}
  for f in concurrent.futures.as_completed(futures):
   result=f.result();d['rows'].append(result);d['rows'].sort(key=lambda r:r['ell']);save();print('completed ell',result['ell'],'seconds',round(result['elapsed_seconds'],2),flush=True)
 d['summary']=summarize(d['rows'],config,cloud,baseline);d['status']='completed_direct_one_sided_contact_audit'
 d['limitations']=['Full production metric through finite ellmax only, with radial one-sided Taylor extrapolation; cutoff and angular convergence must be inspected.','The contact term is diagnosed, not added to the physical source. It should vanish for an exactly matched Lorenz metric.','Bulk Lorenz samples are local checks and do not bound their integral over the whole domain.','Source, radial Green conventions and cloud normalization are shared with the baseline.']
 save();print(json.dumps(d['summary']['cumulative_ell_results'][-1]),flush=True)
if __name__=='__main__':main()
