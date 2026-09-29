import os,sys,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'));OUT=Path(__file__).resolve().parent
from environment_dense_metric import configure_metric_backend
configure_metric_backend(True)
from lorenz_metric import spin1_metric,spin0_metric,nonstatic_metric
from lorenz_ghp import KerrGHP
from lorenz_mode_jet import separated_jet
from li_normalized_cloud import LiNormalizedCloud
os.nice(15);clouds={c:LiNormalizedCloud(ell=c) for c in [1,2]};rows=[];wg=1/(20**1.5+.88)
for r in [2.,5.,15.,25.]:
 t=1.1;g,one=spin1_metric(r,t,20.,.88,1,1,8,return_vector=True,full_current=True);_,zero=spin0_metric(r,t,20.,.88,1,1,8,return_vector=True);_,h=nonstatic_metric(r,t,20.,.88,1,1,order=8);h=np.array([[z.value for z in row] for row in h]);vector=[-1j*(one[i]+zero[i])/wg for i in range(4)]
 for c,cloud in clouds.items():
  R,Rp,Rpp=cloud.radial_state(r)[:,0];S,Sp,A=cloud.angular(t);gc=KerrGHP(r,t,.88,omega=cloud.spectral_omega,m=c,order=8)
  field=separated_jet(gc,0,A+.88**2*cloud.spectral_omega**2-2*.88*c*cloud.spectral_omega,R,Rp,S,Sp,mass_squared=.3**2);grad=[-1j*cloud.omega*field,field.derivative(0),field.derivative(1),1j*c*field]
  for sign in [-1,1]:
   v=vector if sign==1 else [z.conjugate() for z in vector];f=sum(g.inv[i][j]*v[j]*grad[i] for i in range(4) for j in range(4));gf=KerrGHP(r,t,.88,omega=cloud.omega+sign*wg,m=c+sign,order=8)
   terms=[(gf.inv[i][j]*(gf.partial(gf.partial(f,j),i)-sum(gf.gamma[k][i][j]*gf.partial(f,k) for k in range(4)))).value for i in range(4) for j in range(4)];lhs=sum(terms)-.3**2*f.value;rhs=cloud.lorenz_source(r,t,h if sign==1 else h.conjugate());error=abs(lhs-rhs)
   rows.append({'r':r,'theta':t,'cloud':c,'mg':sign,'lhs':[lhs.real,lhs.imag],'rhs':[rhs.real,rhs.imag],'relative_to_source':float(error/max(abs(rhs),1e-300)),'relative_to_operator_terms':float(error/(sum(abs(z) for z in terms)+.3**2*abs(f.value)))})
(OUT/'point_source_identity.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows),flush=True)
