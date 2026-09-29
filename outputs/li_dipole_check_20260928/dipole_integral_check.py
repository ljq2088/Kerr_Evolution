import os,sys,json,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'));OUT=Path(__file__).resolve().parent;BANK=ROOT/'outputs/li_fig9_10_20260928'
from environment_dense_metric import DenseLorenzMetricMode
from li_normalized_cloud import LiNormalizedCloud
from li_field_green import LiFieldGreen
from environment_source import angular_mode,connection
os.nice(15)
metric=DenseLorenzMetricMode(20.,.88,1,1);clouds={c:LiNormalizedCloud(ell=c) for c in [1,2]}
cases=[(1,0,0),(1,2,0),(1,2,2),(2,1,1),(2,3,1),(2,3,3)];data={};greens={}
for key in cases:
 c,l,m=key;d=json.loads((BANK/f'r20_c{c}_l{l}_m{m}.json').read_text())
 if 'origin' in d:d=json.loads(Path(d['origin']).read_text());d['radii']=d['source_radii'];d['weights']=d['source_weights']
 data[key]=d;greens[key]=LiFieldGreen(.88,.3,d['omega'],l,m)
radii=np.array(data[cases[0]]['radii']);weights=np.array(data[cases[0]]['weights']);q=24
assert all(np.array_equal(radii,d['radii']) for d in data.values())
x,w=np.polynomial.legendre.leggauss(q);theta=np.arccos(x);targets={key:angular_mode(theta,key[1],key[2],.88**2*(data[key]['omega']**2-.3**2))[0] for key in cases}
path=OUT/'dipole_source_samples.json';sample=json.loads(path.read_text()) if path.exists() else {'q':q,'cases':cases,'radii':radii.tolist(),'source':[]}
for r in radii[len(sample['source']):]:
 values={key:[] for key in cases}
 for k,t in enumerate(theta):
  h=metric(r,t);inv,gamma=connection(r,t,.88)
  for c,cloud in clouds.items():
   _,H,_=cloud.hessian(r,t);raised=inv@H@inv
   for key in cases:
    if key[0]!=c:continue
    hm=h if key[2]-c==1 else h.conjugate();v=np.einsum('ij,ij->',hm,raised)*(r*r+.88**2*x[k]**2);values[key].append(v)
 row=[2*np.pi*np.dot(w*targets[key],values[key]) for key in cases];sample['source'].append([[z.real,z.imag] for z in row]);tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(sample));tmp.replace(path)
 if len(sample['source'])%8==0:print('radii',len(sample['source']),len(radii),flush=True)
sources=np.array(sample['source']);sources=sources[:,:,0]+1j*sources[:,:,1];identity=json.loads((OUT/'green_identity.json').read_text());rows=[]
for k,key in enumerate(cases):
 g=greens[key];Z=np.sum(weights*g.upsol.sol(radii)[0]*sources[:,k])/g.w0;d=data[key];sj=np.array(d['source']);full=np.sum(weights*g.upsol.sol(radii)[0]*(sj[:,0]+1j*sj[:,1]))/g.w0
 ref=next(v for v in identity['rows'] if (v['cloud'],v['ell'],v['m'])==key and v['quadrature']==40 and v['orbit_offset']==.0005);boundary=complex(*ref['dipole_z'])
 rows.append({'cloud':key[0],'ell':key[1],'m':key[2],'dipole_volume_z':[Z.real,Z.imag],'volume_vs_boundary_amplitude_error':float(abs(Z/boundary-1)),'volume_vs_boundary_flux_difference':float(abs(Z/boundary)**2-1),'dipole_vs_full_amplitude_residual':float(abs(Z/full-1)),'dipole_flux_over_full':float(abs(Z/full)**2)})
(OUT/'dipole_integral_check.json').write_text(json.dumps({'rows':rows,'q':q,'radial_points':len(radii),'limitations':['Finite radial quadrature retained. Same metric but independent direct covariant source and endpoint identity. Does not validate particle matching jumps.']},indent=2));print(json.dumps(rows),flush=True)
