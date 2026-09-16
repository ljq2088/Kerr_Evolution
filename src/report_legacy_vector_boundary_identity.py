"""Independent Green identity for deletion of the dipole vector sector.

Evaluates only the gauge-induced scalar and one-sided orbit/boundary terms,
not the h:Hessian radial quadrature used by report_legacy_vector_ablation.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from environment_source import ThresholdCloud,angular_mode
from environment_radial import RadialGreen
from environment_cloud import mode_flux
from lorenz_metric import spin1_metric
from lorenz_ghp import KerrGHP
from lorenz_mode_jet import separated_jet
from source_provenance import local_dependency_hashes,source_fingerprint,validate_saved_samples
ROOT=Path(__file__).resolve().parents[1]
def enc(z):
 a=np.asarray(z);return np.stack([a.real,a.imag],axis=-1).tolist()
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--baseline',required=True);parser.add_argument('--output',required=True);parser.add_argument('--angular-order',type=int,default=24);parser.add_argument('--jet-order',type=int,default=8);args=parser.parse_args()
 path=Path(args.baseline);d=json.loads(path.read_text());p=d['parameters'];m=p['metric']['m_g'];sm=p['scalar_m'];ell=p['scalar_ell'];r0=p['metric']['orbital_radius']
 if abs(m)!=1 or 'background' in p:raise ValueError('Only threshold Kerr mg+/-1 supported')
 out=Path(args.output)
 if out.exists():raise FileExistsError(out)
 try:validate_saved_samples(d,source_fingerprint());baseline_state='same_current_source_provenance'
 except ValueError as exc:baseline_state='historical_baseline: '+str(exc)
 cloud=ThresholdCloud(alpha=p['alpha'])
 if abs(cloud.a-p['metric']['a'])>1e-14:raise ValueError('cloud mismatch')
 green=RadialGreen(cloud.a,cloud.mu,p['omega'],ell,sm,rmax=p['green_outer_radius'],offset=p['green_horizon_offset'],rtol=1e-11,infinity_method=p.get('infinity_method','series'))
 x,w=np.polynomial.legendre.leggauss(args.angular_order);theta=np.arccos(x);angular=angular_mode(theta,ell,sm,cloud.a**2*(p['omega']**2-cloud.mu**2))[0];samples=[]
 def state(r):
  R,Rp=cloud.radial(r)[:,0];values=[]
  for t in theta:
   g,xi=spin1_metric(r,t,r0,cloud.a,1,1,args.jet_order,return_vector=True,full_current=True)
   minus_lie=[-1j*v/g.omega for v in xi]
   if m<0:minus_lie=[v.conjugate() for v in minus_lie]
   gc=KerrGHP(r,t,cloud.a,omega=cloud.omega,m=1,order=args.jet_order)
   S,Sp,A=angular_mode(t,1,1,cloud.c2);lam=A+cloud.a**2*cloud.omega**2-2*cloud.a*cloud.omega
   field=separated_jet(gc,0,lam,R,Rp,S,Sp,mass_squared=cloud.mu**2)
   f=sum(g.inv[i][j]*minus_lie[j]*gc.partial(field,i) for i in range(4) for j in range(4))
   values.append([f.value,f.derivative(0).value,f.derivative(0).derivative(0).value])
  projected=2*np.pi*(w*angular)@np.asarray(values);samples.append(dict(r=float(r),state=enc(projected)));return projected
 def boundary(r,st,kind):
  v,dv=(green.upsol if kind=='horizon' else green.insol).sol(r)
  return (r*r-2*r+cloud.a**2)*(v*st[1]-dv*st[0])/green.w0
 inner,outer=p['source_panels'][0],p['source_panels'][-1];endstates=[state(inner),state(outer)]
 rows=[];fullH=complex(*d['z_h']);fullI=complex(*d['z_inf'])
 for eps in [1e-3,5e-4]:
  limits=[]
  for sign in [-1,1]:
   st=state(r0+sign*eps);dr=-sign*eps;limits.append([st[0]+dr*st[1]+.5*dr*dr*st[2],st[1]+dr*st[2]])
  result={}
  for kind in ['horizon','infinity']:
   endpoint=boundary(outer,endstates[1],kind)-boundary(inner,endstates[0],kind)
   orbit=boundary(r0,limits[0],kind)-boundary(r0,limits[1],kind)
   result[kind]=dict(endpoint=enc(endpoint),orbit=enc(orbit),Zvector=enc(endpoint+orbit))
  zvh=complex(*result['horizon']['Zvector']);zvi=complex(*result['infinity']['Zvector']) if green.propagating else 0j
  legacyH=fullH-zvh;legacyI=fullI-zvi
  rows.append(dict(orbit_extrapolation_epsilon=eps,boundaries=result,z_h_legacy=enc(legacyH),z_inf_legacy=enc(legacyI),flux_legacy=mode_flux(p['omega'],sm,cloud.omega,cloud.m,cloud.mu,cloud.a,legacyI,legacyH),vector_horizon_amplitude_over_full=enc(zvh/fullH)))
 record=dict(status='completed_vacuum_vector_boundary_identity',baseline=str(path),baseline_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),baseline_provenance=baseline_state,parameters=p,angular_order=args.angular_order,jet_order=args.jet_order,implementation_sha256=local_dependency_hashes(ROOT/'src',['report_legacy_vector_boundary_identity']),rows=rows,samples=samples,baseline_z_h=d['z_h'],baseline_flux=d['flux'],limitations=['Pure-gauge scalar identity on each vacuum side, including one-sided orbit terms; these terms are NOT new physical delta sources.','Tests a numerical deletion hypothesis, not an admissible alternative sourced metric or attribution of the2025 figure.','Full baseline retained as recorded; historical provenance is explicitly marked.'])
 out.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n');print(json.dumps(rows,indent=2),flush=True)
if __name__=='__main__':main()
