"""Summarize orbit contact convergence without recomputing metric data."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from environment_source import ThresholdCloud,angular_mode
from lorenz_ghp import KerrGHP

def enc(z):
 a=np.asarray(z,complex);return np.stack([a.real,a.imag],axis=-1).tolist()
def dec(z):
 a=np.asarray(z);return a[...,0]+1j*a[...,1]
def main():
 parser=argparse.ArgumentParser();parser.add_argument('inputs',nargs='+');parser.add_argument('--output',required=True);args=parser.parse_args();results=[]
 for name in args.inputs:
  path=Path(name);d=json.loads(path.read_text())
  if not d['status'].startswith('completed'):raise ValueError('Unfinished input')
  c=d['metadata']['config'];cloud=ThresholdCloud(alpha=c['alpha']);r0=c['r0'];a=cloud.a;q=c['quadrature'];eps=c['epsilon'];x,w=np.polynomial.legendre.leggauss(q);tt=np.arccos(x);omega=d['summary']['omega'];S=angular_mode(tt,0,0,a*a*(omega*omega-cloud.mu*cloud.mu))[0];Sb,Sbp,_=angular_mode(tt,1,1,cloud.c2)
  cumulative=[];constraint=np.zeros((2,q,4),complex)
  for row in d['rows']:
   constraint+=dec(row['lorenz_constraint_covariant_at_samples']);bulk=[]
   for index,sign in enumerate((-1,1)):
    r=r0+sign*eps;R,Rp=cloud.radial(r)[:,0];grad=np.array([-1j*cloud.omega*R*Sb,Rp*Sb,R*Sbp,1j*R*Sb]).T
    inv=np.asarray([[[z.value for z in rr] for rr in KerrGHP(r,t,a,order=2).inv] for t in tt]);integrand=(r*r+a*a*x*x)*np.einsum('kab,kb,ka->k',inv,constraint[index],grad)
    bulk.append(2*np.pi*np.dot(w*S,integrand))
   cumulative.append(dict(ellmax=row['ell'],bulk_C_grad_phi_scalar00_projection=enc(bulk),max_covariant_constraint_component=float(abs(constraint).max())))
  results.append(dict(input=str(path),input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),config=c,endpoint=d['summary']['cumulative_ell_results'][-1],contact_by_ell=d['summary']['cumulative_ell_results'],bulk_by_ell=cumulative,source_provenance=d['metadata']['source_provenance']['sha256']))
 record=dict(status='contact_projection_convergence_summary_not_global_Lorenz_proof',results=results,interpretation='deltaZH is the Green response to the directly measured jump [bar h^{rb}] grad_b Phi0; local bulk samples are not a radial integral bound',implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 Path(args.output).write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
 for r in results:print(r['input'],r['endpoint']['relative_delta_z_h'],r['endpoint']['signed_flux_relative_change'],r['bulk_by_ell'][-1])
if __name__=='__main__':main()
