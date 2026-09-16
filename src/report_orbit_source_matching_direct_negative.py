"""Direct m_g=-1 control for the orbit-contact diagnostic.

Reuses only persistence/projection from the conjugate-route driver. This worker
calls the negative-frequency metric directly and never conjugates its output.
"""
import sys,time,json,hashlib
from pathlib import Path
import numpy as np
import report_orbit_source_matching_contact as common
from lorenz_metric import nonstatic_metric
from lorenz_tensor import lorenz_constraint

def direct_work(task):
 ell,c=task;start=time.monotonic();r0=c['r0'];a=c['a'];eps=c['epsilon'];q=c['quadrature'];order=c['jet_order']
 theta=np.arccos(np.polynomial.legendre.leggauss(q)[0]);limits=[];linear=[];bulk=[]
 for sign in (-1,1):
  r=r0+sign*eps;dr=-sign*eps;side=[];sidelinear=[];sidebulk=[]
  for t in theta:
   g,h=nonstatic_metric(r,t,r0,a,ell,-1,order=order)
   v=np.asarray([[z.value for z in row] for row in h]);v1=np.asarray([[z.derivative(0).value for z in row] for row in h]);v2=np.asarray([[z.derivative(0).derivative(0).value for z in row] for row in h])
   side.append(v+dr*v1+dr*dr*v2/2);sidelinear.append(v+dr*v1)
   sidebulk.append(np.asarray([z.value for z in lorenz_constraint(g,h)]))
  limits.append(side);linear.append(sidelinear);bulk.append(sidebulk)
 return dict(ell=ell,metric_limits=common.enc(limits),metric_limits_linear=common.enc(linear),lorenz_constraint_covariant_at_samples=common.enc(bulk),elapsed_seconds=time.monotonic()-start)

if __name__=='__main__':
 output=Path(sys.argv[sys.argv.index('--output')+1])
 if output.exists():raise FileExistsError('Direct-negative control requires a new output; never resume a conjugate-route checkpoint')
 common.work=direct_work
 common.main()
 d=json.loads(output.read_text());d['route']='direct_nonstatic_metric_m_minus_one_no_conjugation';d['driver_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();output.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
