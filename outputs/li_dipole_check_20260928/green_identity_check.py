import os,sys,json,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'));OUT=Path(__file__).resolve().parent;BANK=ROOT/'outputs/li_fig9_10_20260928'
from environment_dense_metric import configure_metric_backend
configure_metric_backend(True)
from lorenz_metric import spin1_metric,spin0_metric
from lorenz_mode_jet import separated_jet
from lorenz_ghp import KerrGHP
from li_normalized_cloud import LiNormalizedCloud
from li_field_green import LiFieldGreen
from environment_source import angular_mode
os.nice(15)
def enc(v):
 a=np.asarray(v,complex);return np.stack([a.real,a.imag],axis=-1).tolist()
clouds={c:LiNormalizedCloud(ell=c) for c in [1,2]}
cases=[(1,0,0),(1,2,0),(1,2,2),(2,1,1),(2,3,1),(2,3,3)]
data={};greens={};fullz={}
for key in cases:
 c,l,m=key;d=json.loads((BANK/f'r20_c{c}_l{l}_m{m}.json').read_text())
 if 'origin' in d:
  d=json.loads(Path(d['origin']).read_text());d['radii']=d['source_radii'];d['weights']=d['source_weights']
 data[key]=d;g=LiFieldGreen(.88,.3,d['omega'],l,m);greens[key]=g
 rr=np.array(d['radii']);w=np.array(d['weights']);z=np.array(d['source']);J=z[:,0]+1j*z[:,1];fullz[key]=sum(w*g.upsol.sol(rr)[0]*J)/g.w0
wg=1/(20**1.5+.88);rows=[]
for q in [24,40]:
 x,weights=np.polynomial.legendre.leggauss(q);theta=np.arccos(x);targets={key:angular_mode(theta,key[1],key[2],.88**2*(data[key]['omega']**2-.3**2))[0] for key in cases}
 def state(r):
  values={key:[] for key in cases}
  for t in theta:
   g,one=spin1_metric(r,t,20.,.88,1,1,8,return_vector=True,full_current=True);_,zero=spin0_metric(r,t,20.,.88,1,1,8,return_vector=True)
   vector=[-1j*(one[i]+zero[i])/wg for i in range(4)]
   for c,cloud in clouds.items():
    R,Rp,Rpp=cloud.radial_state(r)[:,0];S,Sp,A=cloud.angular(t);gc=KerrGHP(r,t,.88,omega=cloud.spectral_omega,m=c,order=8)
    field=separated_jet(gc,0,A+.88**2*cloud.spectral_omega**2-2*.88*c*cloud.spectral_omega,R,Rp,S,Sp,mass_squared=.3**2)
    grad=[-1j*cloud.omega*field,field.derivative(0),field.derivative(1),1j*c*field]
    for sign in [-1,1]:
     v=vector if sign==1 else [z.conjugate() for z in vector]
     f=sum(g.inv[i][j]*v[j]*grad[i] for i in range(4) for j in range(4));st=[f.value,f.derivative(0).value,f.derivative(0).derivative(0).value]
     for key in cases:
      if key[0]==c and key[2]-c==sign:values[key].append(st)
  return {key:2*np.pi*(weights*targets[key])@np.array(values[key]) for key in cases}
 def B(key,r,st):
  g=greens[key];u,du=g.upsol.sol(r);return (r*r-2*r+.88**2)*(u*st[1]-du*st[0])/g.w0
 inner=float(data[cases[0]]['radii'][0]);outer=float(data[cases[0]]['radii'][-1])
 # True integration-panel endpoints, not the first/last Gauss nodes.
 inner=clouds[1].rp+5e-4;outer=320.
 states={}
 for r in [inner,outer,19.999,20.001,19.9995,20.0005]:
  states[r]=state(r);print('sampled',q,r,flush=True)
 for key in cases:
  for eps in [.001,.0005]:
   limits=[]
   for sign in [-1,1]:
    st=states[20+sign*eps][key];dr=-sign*eps;limits.append([st[0]+dr*st[1]+.5*dr*dr*st[2],st[1]+dr*st[2]])
   endpoint=B(key,outer,states[outer][key])-B(key,inner,states[inner][key]);jump=B(key,20,limits[0])-B(key,20,limits[1]);z=endpoint+jump;ratio=z/fullz[key]
   rows.append(dict(cloud=key[0],ell=key[1],m=key[2],quadrature=q,orbit_offset=eps,dipole_z=enc(z),endpoint=enc(endpoint),orbit_jump=enc(jump),full_z=enc(fullz[key]),dipole_over_full=enc(ratio),residual_amplitude_fraction=float(abs(1-ratio)),dipole_flux_over_full=float(abs(ratio)**2)))
 (OUT/'green_identity.json').write_text(json.dumps({'status':'boundary_identity_dipole_attribution_not_jump_validation','rows':rows,'conventions':'h_ab=2 nabla_(a v_b), induced scalar=v^a partial_a Phi; complex spatial cloud with real temporal frequency. Orbit jump retained.','limitations':['Full baseline is finite quadrature source integral; identity uses dipole only.','Vacuum identity does not prove particle matching amplitudes correct.','No author data used or normalization fitted.']},indent=2))
 print('RESULT',q,json.dumps(rows[-12:]),flush=True)
