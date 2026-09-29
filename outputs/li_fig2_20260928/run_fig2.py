"""Li v2 Figs 9/10: unadjusted effective scalar fluxes, resumable finite runs."""
import os,sys,json,hashlib,subprocess,traceback,argparse
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import pymupdf as fitz
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from source_provenance import local_dependency_hashes
BANK=ROOT/'docs/li_alignment/rp20_cloud11_j18_L20_q40_nr12_h64_cache512'
ORBIT_PLAN=[20.,18.3,4.1,10.1,30.1,41.1,42.1,50.1,18.1,18.5,40.1]
F10={1:{'infinity':[(3,3),(4,4),(5,5),(2,2)],'horizon':[(0,0),(2,2),(2,0),(1,-1)]},2:{'infinity':[(4,4),(5,5),(5,3),(3,3)],'horizon':[(1,1),(3,3),(3,1),(0,0)]}}
F9=[(l,m) for l in range(2,11) for m in range(2 if l%2==0 else 3,l+1,2)]
def stamp():return datetime.now(timezone.utc).isoformat()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
 p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix('.tmp');t.write_text(json.dumps(d,indent=2,allow_nan=False));t.replace(p)
def encode(z):
 z=np.asarray(z,complex);return np.stack([z.real,z.imag],axis=-1).tolist()
def fingerprint():
 d=local_dependency_hashes(ROOT/'src',['li_normalized_cloud','li_separated_source','li_field_green','report_li_aligned_flux','environment_metric_sampling'])
 d['li_source_terms.json']=sha(ROOT/'src/li_source_terms.json');return d
def check_hashes(d):
 for k,v in d.items():
  p=ROOT/'src'/(k if k.endswith('.json') else k+'.py')
  if not p.exists() or sha(p)!=v:raise ValueError(f'Implementation changed: {k}')
def modepath(r,c,l,m):return OUT/f'r{r:g}_c{c}_l{l}_m{m}.json'
EXPECTED=[(l,m) for l in range(6) for m in range(-l,l+1) if (l+m)%2==0 and m!=1]
ORBIT_PLAN=[20.,18.3,4.1,10.1,30.1,26.1,26.7,27.1,40.1,41.1,41.5,41.6,41.7,41.8,42.1,50.1]
SHARED=ROOT/'outputs/li_fig9_10_20260928'
def targets(r,mg):return [(1,l,m) for l,m in EXPECTED if abs(m-1)==mg]
def group(r,mg,workers):
 from li_normalized_cloud import LiNormalizedCloud
 from li_separated_source import LiSeparatedSource
 from li_field_green import LiFieldGreen
 from environment_source import angular_mode
 from environment_cloud import mode_flux
 from report_environment_forced_mode import source_grid
 from report_li_aligned_flux import spherical_coefficients
 from environment_dense_metric import DenseLorenzMetricMode
 from environment_metric_sampling import precompute_metric
 from types import SimpleNamespace
 seed_cache(r,mg)
 os.nice(10);fp=fingerprint();todo=[t for t in targets(r,mg) if not modepath(r,*t).exists()]
 if not todo:return
 clouds={c:LiNormalizedCloud(ell=c) for c in sorted({v[0] for v in todo})}
 cloud=next(iter(clouds.values()));grid=SimpleNamespace(orbital_radius=r,source_inner_offset=5e-4,source_outer_radius=320.,radial_order=12 if r==20. else 24,horizon_order=64,horizon_log=True)
 panels,radii,weights=source_grid(grid,cloud)
 # Refine rapid oscillations at small orbital radii without changing existing rp20 caches.
 if r!=20.:
  wg=mg/(r**1.5+.88);k=np.sqrt(max((.3+wg)**2-.3**2,0));width=min(40.,12*np.pi/max(wg+k,.01));refined=[panels[0],panels[1]]
  for lo,hi in zip(panels[1:-1],panels[2:]):refined.extend(np.linspace(lo,hi,int(np.ceil((hi-lo)/width))+1)[1:])
  panels=np.array(refined);rr=[];ww=[]
  for i,(lo,hi) in enumerate(zip(panels[:-1],panels[1:])):
   x,w=np.polynomial.legendre.leggauss(64 if i==0 else 24)
   if i==0:
    a,b=np.log(lo-cloud.rp),np.log(hi-cloud.rp);dist=np.exp((a+b)/2+(b-a)/2*x);rr.extend(cloud.rp+dist);ww.extend((b-a)/2*w*dist)
   else:rr.extend((lo+hi)/2+(hi-lo)/2*x);ww.extend((hi-lo)/2*w)
  radii=np.array(rr);weights=np.array(ww)
 x,wt=np.polynomial.legendre.leggauss(40);theta=np.arccos(x);banks={};audits={};needed=sorted({m-c for c,l,m in todo})
 for signed in needed:
  old=BANK/f'metric_spherical_mg{signed}.npz'
  if r==20. and old.exists():
   with np.load(old,allow_pickle=False) as b:
    meta=json.loads(str(b['metadata']));check_hashes(meta['implementation_sha256']);assert np.array_equal(b['radii'],radii);banks[signed]=b['coefficients'].copy()
   audits[signed]={'reused':str(old),'sha256':sha(old)}
 if len(banks)!=len(needed):
  metric,audit=precompute_metric(DenseLorenzMetricMode(r,.88,mg,20),radii,theta,OUT/'metric_cache',workers=workers)
  for signed in needed:
   if signed in banks:continue
   banks[signed]=np.array([spherical_coefficients(metric,z,theta,wt,.88,signed,18) for z in radii]);audits[signed]=audit
   bp=OUT/f'metric_r{r:g}_mg{signed}.npz';np.savez_compressed(bp,radii=radii,coefficients=banks[signed],metadata=json.dumps({'implementation_sha256':fp,'audit':audit}));audits[signed]=dict(audit,bank_sha256=sha(bp))
 for c,l,m in todo:
  cloud=clouds[c];signed=m-c;omega=cloud.omega+signed/(r**1.5+.88)
  projector=LiSeparatedSource(a=.88,omega_c=cloud.omega,m_c=c,metric_m=signed,spherical_lmax=18,cloud_angular=cloud.angular_state,target_angular=lambda t:angular_mode(t,l,m,.88**2*(omega**2-.3**2))[0],quadrature=64,pmax=12,dps=64)
  checkpoint=OUT/f'checkpoint_r{r:g}_c{c}_l{l}_m{m}.json';J=[]
  if checkpoint.exists():
   saved=json.loads(checkpoint.read_text());assert saved['implementation_sha256']==fp and saved['radii']==radii.tolist();J=[complex(*z) for z in saved['source']]
  for z,coef in zip(radii[len(J):],banks[signed][len(J):]):
   J.append(projector.project(z,cloud.radial_state(z)[:,0],coef));save(checkpoint,{'implementation_sha256':fp,'radii':radii.tolist(),'source':encode(J)})
  J=np.asarray(J);green=LiFieldGreen(.88,.3,omega,l,m);u=green.insol.sol(radii)[0];v=green.upsol.sol(radii)[0];zi=np.sum(weights*u*J)/green.w0;zh=np.sum(weights*v*J)/green.w0
  spread=float(np.max(abs(green.wronskian(radii)/green.w0-1)))
  if not np.isfinite([zi,zh]).all() or spread>1e-5:raise ValueError(f'Radial closure failed {spread}')
  flux=mode_flux(omega,m,cloud.omega,c,.3,.88,zi,zh);om=1/(r**1.5+.88)
  for b in flux:
   if abs(flux[b]['orbital_energy']-om*flux[b]['orbital_angular_momentum'])>1e-12*max(abs(flux[b]['orbital_energy']),1e-30):raise ValueError('E=Omega L failed')
  cutoffs=[]
  for end in (80.,160.,320.):
   sel=radii<end;fz=mode_flux(omega,m,cloud.omega,c,.3,.88,np.sum(weights[sel]*u[sel]*J[sel])/green.w0,np.sum(weights[sel]*v[sel]*J[sel])/green.w0);cutoffs.append({'source_outer':end,'flux':{b:fz[b]['orbital_energy']/.3**6 for b in fz}})
  row={'status':'new_finite_resolution_mode','r0':r,'cloud':c,'ell':l,'m':m,'omega':omega,'cloud_provenance':cloud.provenance,'implementation_sha256':fp,'metric':audits[signed],'radii':radii.tolist(),'weights':weights.tolist(),'source':encode(J),'flux':{b:flux[b]['orbital_energy']/.3**6 for b in flux},'flux_details_unit_mass':flux,'source_outer_cutoff_sequence':cutoffs,'infinity_cancellation_condition':float(np.sum(abs(weights*u*J))/max(abs(np.sum(weights*u*J)),1e-300)),'wronskian_relative_spread':spread,'boundary':green.boundary_audit,'green_outer':green.rmax,'completed_utc':stamp()}
  save(modepath(r,c,l,m),row);print('COMPLETED',r,c,l,m,row['flux'],flush=True);render()


import shutil,time

def seed_cache(r,mg):
 # Snapshot committed files only. Never write into the running Fig10 cache.
 source=ROOT/'outputs/metric_cache';dest=OUT/'metric_cache'
 copied=0
 for folder in source.iterdir():
  files=list(folder.glob('*.npz'))
  if not files:continue
  with np.load(files[0],allow_pickle=False) as v:meta=json.loads(str(v['metadata']))
  m=meta['metric']
  if m.get('orbital_radius')!=r or m.get('m_g')!=mg:continue
  target=dest/folder.name;target.mkdir(parents=True,exist_ok=True)
  for f in files:
   if not (target/f.name).exists():shutil.copy2(f,target/f.name);copied+=1
 print('Seeded committed metric radii',r,mg,copied,flush=True)

def ingest():
 for p in SHARED.glob('r*_c1_l*_m*.json'):
  d=json.loads(p.read_text())
  if (d['ell'],d['m']) not in EXPECTED:continue
  check_hashes(d['implementation_sha256']);dst=modepath(d['r0'],1,d['ell'],d['m'])
  if not dst.exists():save(dst,dict(d,reused_file=str(p),reused_sha256=sha(p)))

def reference():
 pdf=ROOT/'outputs/paper_original_reference/li_2507_02045v2/source/total_flux_11.pdf'
 calibration=ROOT/'docs/environment_reproduction/li_2025_figure2_readout_20260916.json'
 expected=json.loads(calibration.read_text());assert sha(pdf)==expected['pdf_sha256']
 doc=fitz.open(pdf);page=doc[0];ds=page.get_drawings()
 nx=np.polyfit([61.12561798095703,98.98464965820312,136.8436737060547,174.7027130126953,212.56173706054688],[10,20,30,40,50],1)
 ny=np.polyfit([148.2509307861328,110.11843872070312,71.98595428466797,33.85345458984375],[-4,-3,-2,-1],1)
 curves={}
 for idx,name in [(41,'li_horizon'),(42,'li_infinity'),(43,'dyson_horizon'),(44,'dyson_infinity')]:
  items=ds[idx]['items'];assert all(i[0]=='l' for i in items)
  pts=np.array([tuple(items[0][1])]+[tuple(i[2]) for i in items]);rr=np.polyval(nx,pts[:,0]);vv=10**np.polyval(ny,pts[:,1]);assert np.all(np.diff(rr)>0)
  curves[name]={'r':rr.tolist(),'magnitude':vv.tolist()}
 save(OUT/'reference.json',{'pdf_sha256':sha(pdf),'curves':curves,'type':'Original published vector curves, not author numerical arrays'})
 page.get_pixmap(matrix=fitz.Matrix(3,3)).save(OUT/'paper_fig2.png')

def render():
 ref=json.loads((OUT/'reference.json').read_text());curves=ref['curves'];rows=[];coverage=[]
 for r in ORBIT_PLAN:
  modes=[]
  for l,m in EXPECTED:
   p=modepath(r,1,l,m)
   if p.exists():modes.append(json.loads(p.read_text()))
  coverage.append({'r0':r,'completed':len(modes),'expected':len(EXPECTED)})
  if len(modes)!=len(EXPECTED):continue
  for b in ('infinity','horizon'):
   v=sum(d['flux'][b] for d in modes);c=curves['li_'+b];paper=float(10**np.interp(r,c['r'],np.log10(c['magnitude'])))
   rows.append({'r0':r,'boundary':b,'computed_signed_Li_units':v,'paper_magnitude':paper,'difference_percent':100*(abs(v)/paper-1)})
 save(OUT/'comparison.json',{'rows':rows,'coverage':coverage,'normalization':'Effective orbital energy flux per epsilon^2 zeta^2. No fitted multipliers.','limitations':['Only all-18-mode totals plotted. Finite resolution, convergence pending.','Dashed curves are published Dyson comparison curves, not a new independent Dyson run.','Paper values digitized; interpolation near thresholds needs caution.']})
 fig,(ax,err)=plt.subplots(2,1,figsize=(10,7.2),sharex=True,layout='constrained',gridspec_kw={'height_ratios':[3,1]})
 for b,color,label in [('infinity','#1976b8','Infinity'),('horizon','#db7a18','Minus horizon')]:
  for name,style in [('li','-'),('dyson','--')]:
   c=curves[name+'_'+b];ax.semilogy(c['r'],c['magnitude'],style,color=color,lw=1.6,alpha=1 if name=='li' else .65,label=f'{label}: '+('Li' if name=='li' else 'Dyson as plotted by Li'))
  data=[v for v in rows if v['boundary']==b]
  ax.semilogy([v['r0'] for v in data],[abs(v['computed_signed_Li_units']) for v in data],'o',mfc='white',mec=color,mew=1.8,ms=8,label=f'{label}: our complete sums')
  err.plot([v['r0'] for v in data],[v['difference_percent'] for v in data],'o',color=color,ms=6)
 ax.axvline(41.66,color='#bb3344',lw=1,alpha=.7);ax.set(ylabel=r'$|\dot E^s|/(\epsilon^2\zeta^2)$',xlim=(4.1,50.1));ax.legend(fontsize=8,ncol=2);ax.grid(alpha=.2)
 err.axhline(0,color='gray',lw=.8);err.set(xlabel=r'$r_0/M$',ylabel='Our vs Li [%]');err.grid(alpha=.2)
 fig.suptitle('Li Fig.2 comparison | a/M=0.88, M mu=0.3, cloud |211>, ell <= 5\nOnly completed numerical points shown; no fitted normalization',fontsize=12)
 fig.savefig(OUT/'figure2_comparison.png',dpi=180);plt.close(fig)

def run():
 ingest();render();fp=fingerprint()
 previous=json.loads((OUT/'execution.json').read_text()) if (OUT/'execution.json').exists() else {}
 if previous and previous['implementation_sha256']!=fp:raise RuntimeError('Numerical implementation changed')
 for f in OUT.glob('r*_c1_l*_m*.json'):check_hashes(json.loads(f.read_text())['implementation_sha256'])
 state={'status':'running','pid':os.getpid(),'started_utc':stamp(),'implementation_sha256':fp,'orbit_plan':ORBIT_PLAN,'workers':4,'completed_groups':[],'failures':[]}
 save(OUT/'execution.json',state)
 for r in ORBIT_PLAN:
  for mg in [1,5,6,2,3,4]:
   ingest()
   if all(modepath(r,*t).exists() for t in targets(r,mg)):continue
   # Avoid recalculating the Fig10 group's in-flight radii; its durable parent continues independently.
   while True:
    other=json.loads((SHARED/'execution.json').read_text());pid=other.get('pid')
    busy=other.get('status')=='running' and pid and Path(f'/proc/{pid}').exists() and other.get('active_orbit')==r and other.get('active_metric_abs_m')==mg
    if not busy:break
    state.update(status='waiting_for_shared_metric',active_orbit=r,active_metric_abs_m=mg,updated_utc=stamp());save(OUT/'execution.json',state);time.sleep(60)
   ingest();state.update(status='running',active_orbit=r,active_metric_abs_m=mg,updated_utc=stamp());save(OUT/'execution.json',state)
   assert fingerprint()==fp,'Numerical sources changed'
   with (OUT/f'group_r{r:g}_mg{mg}.log').open('a') as log:
    proc=subprocess.run([sys.executable,'-u',__file__,'--group',str(r),str(mg)],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
   if proc.returncode:state['failures'].append({'r0':r,'mg':mg,'returncode':proc.returncode})
   else:state['completed_groups'].append([r,mg])
   save(OUT/'execution.json',state);render()
 state.update(status='partial_with_failures' if state['failures'] else 'completed_finite_scan_convergence_pending',completed_utc=stamp());save(OUT/'execution.json',state)

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--run',action='store_true');ap.add_argument('--group',nargs=2,type=float);args=ap.parse_args()
 if args.group:group(args.group[0],int(args.group[1]),4)
 else:
  reference();ingest();render()
  if args.run:
   try:run()
   except Exception:
    p=OUT/'execution.json';d=json.loads(p.read_text());d.update(status='failed_resumable',traceback=traceback.format_exc());save(p,d);raise
