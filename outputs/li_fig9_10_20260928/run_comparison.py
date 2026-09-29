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
def targets(r,mg):
 result=[]
 for c in (1,2):
  cases=set(sum(F10[c].values(),[]))
  if c==1 and r==20.:cases.update(F9)
  result.extend((c,l,m) for l,m in sorted(cases) if abs(m-c)==mg)
 return result

def reference():
 pdf=ROOT/'outputs/paper_original_reference/li_2507_02045v2/source/mode_by_mode_flux.pdf';p=fitz.open(pdf)[0];d=p.get_drawings()
 xt=[d[i]['rect'].x0 for i in [26,29,32,35,38]];xm=np.polyfit(xt,[2,4,6,8,10],1)
 ym=np.polyfit([d[i]['rect'].y0 for i in [41,44,47,50,53,56]],[-13,-11,-9,-7,-5,-3],1)
 rows=[]
 for l in range(2,11):
  items=[z for z in d[2:26] if round(np.polyval(xm,(z['rect'].x0+z['rect'].x1)/2))==l]
  items.sort(key=lambda z:z['rect'].y0)
  for rank,z in enumerate(items):
   m=l-2*rank;y=(z['rect'].y0+z['rect'].y1)/2
   rows.append({'ell':l,'m':m,'Li':float(10**np.polyval(ym,y))})
 p10=ROOT/'outputs/paper_original_reference/li_2507_02045v2/source/flux_dominant_modes.pdf';d=fitz.open(p10)[0].get_drawings();curves=[]
 for c,first,rectidx,ymin,ymax in [(1,38,1,-4,2),(2,111,60,-3,6)]:
  # Plot bounds and grid lines, using original vector coordinates rather than pixels.
  left=50.801666259765625 if c==1 else 287.5159606933594
  grid=[v['rect'] for v in d if v['color'] and abs(v['color'][0]-.6901198)<1e-4 and v['rect'].width<1e-3 and v['rect'].height>100 and (v['rect'].x0<260 if c==1 else v['rect'].x0>280)]
  xx=sorted(set(v.x0 for v in grid));assert len(xx)==5
  xm=np.polyfit(xx,[10,20,30,40,50],1)
  horizontals=[v['rect'] for v in d if v['color'] and abs(v['color'][0]-.6901198)<1e-4 and v['rect'].height<1e-3 and v['rect'].width>150 and (v['rect'].x0<260 if c==1 else v['rect'].x0>280)]
  yy=sorted(set(v.y0 for v in horizontals));ym=np.polyfit([yy[0],yy[-1]],[ymax,ymin],1)
  for boundary,offset in [('horizon',0),('infinity',1)]:
   for j,(l,m) in enumerate(F10[c][boundary]):
    items=d[first+offset+2*j]['items'];assert all(z[0]=='l' for z in items)
    a=np.array([tuple(items[0][1])]+[tuple(z[2]) for z in items]);r=np.polyval(xm,a[:,0]);v=10**np.polyval(ym,a[:,1]);assert np.all(np.diff(r)>0)
    curves.append({'cloud':c,'ell':l,'m':m,'boundary':boundary,'r':r.tolist(),'Li_magnitude':v.tolist()})
 result={'fig9':rows,'fig10':curves,'sha256':{pdf.name:sha(pdf),p10.name:sha(p10)},'limitations':['Vector plot readout, not author raw arrays; no rigorous digitization error.','Fig9 m labels assigned by parity and descending m from top to bottom; ell10,m2 is not visible.']}
 save(OUT/'reference.json',result);return result

def ingest():
 for p in BANK.glob('mode_l*_m*.json'):
  d=json.loads(p.read_text());check_hashes(d['implementation_sha256']);c=d['config']
  assert c['a']==.88 and c['r0']==20. and c['cloud_ell']==1
  dest=modepath(20.,1,d['ell'],d['m'])
  if dest.exists():continue
  save(dest,{'status':'reused_finite_resolution_result','r0':20.,'cloud':1,'ell':d['ell'],'m':d['m'],'flux':{b:d['flux'][b]['orbital_energy']/.3**6 for b in ('infinity','horizon')},'origin':str(p),'origin_sha256':sha(p),'implementation_sha256':d['implementation_sha256'],'wronskian_relative_spread':d['wronskian_relative_spread'],'boundary':d['infinity_boundary'],'source_outer_cutoff_sequence':d['source_outer_cutoff_sequence']})

def render():
 refs=json.loads((OUT/'reference.json').read_text());records=[json.loads(p.read_text()) for p in OUT.glob('r*_c*_l*_m*.json')];lookup={(d['r0'],d['cloud'],d['ell'],d['m']):d for d in records};comp9=[];comp10=[]
 fig,ax=plt.subplots(2,1,figsize=(8,7),sharex=True,layout='constrained',gridspec_kw={'height_ratios':[3,1]})
 for ref in refs['fig9']:
  l,m=ref['ell'],ref['m'];ax[0].semilogy(l,ref['Li'],'x',color='C1')
  row=lookup.get((20.,1,l,m))
  if row:
   v=row['flux']['infinity'];ratio=v/ref['Li'];ax[0].semilogy(l,v,'o',mfc='none',color='C0');ax[1].plot(l,100*(ratio-1),'o',ms=4,color='C0');comp9.append(dict(ell=l,m=m,local=v,Li=ref['Li'],relative_difference_percent=100*(ratio-1)))
 ax[0].plot([],[],'x',color='C1',label='Li Fig.9');ax[0].plot([],[],'o',mfc='none',color='C0',label='Local: no fitted rescaling');ax[0].legend();ax[0].set(ylabel='Infinity effective scalar flux',title=f'Li Fig.9: {len(comp9)}/{len(refs["fig9"])} visible points compared');ax[0].grid(alpha=.2)
 ax[1].axhline(0,color='k',lw=.8);ax[1].set(xlabel='ell',ylabel='Difference [%]');ax[1].grid(alpha=.2);fig.savefig(OUT/'figure9_comparison.png',dpi=170);plt.close(fig)
 fig,axs=plt.subplots(2,2,figsize=(13,8),sharex='col',layout='constrained',gridspec_kw={'height_ratios':[3,1]})
 for c in (1,2):
  axes=axs[:,c-1]
  for ref in refs['fig10']:
   if ref['cloud']!=c:continue
   l,m,b=ref['ell'],ref['m'],ref['boundary'];j=F10[c][b].index((l,m));color='C0' if b=='infinity' else 'C1';ls=['-','--','-.',':'][j]
   rr=np.array(ref['r']);vv=np.array(ref['Li_magnitude']);axes[0].semilogy(rr,vv,ls,color=color,lw=1.2,label=f'{b} ({l},{m})')
   for r in ORBIT_PLAN:
    row=lookup.get((r,c,l,m))
    if row is None or not rr.min()<=r<=rr.max():continue
    v=row['flux'][b];reference=float(10**np.interp(r,rr,np.log10(vv)));ratio=abs(v)/reference
    if v!=0:axes[0].semilogy(r,abs(v),'o',mfc='white',ms=3,color=color)
    axes[1].plot(r,100*(ratio-1),'o',ms=3,color=color)
    comp10.append(dict(cloud=c,r0=r,ell=l,m=m,boundary=b,local_signed=v,Li_magnitude=reference,relative_magnitude_difference_percent=100*(ratio-1)))
  axes[0].set(title=f'Cloud ell_c=m_c={c}',ylabel='|Effective scalar energy flux|',xlim=(4.1,50.1));axes[0].legend(fontsize=7,ncol=2);axes[0].grid(alpha=.2);axes[1].axhline(0,color='k',lw=.8);axes[1].set(xlabel='r0/M',ylabel='Difference [%]');axes[1].grid(alpha=.2)
 fig.suptitle('Li Fig.10: lines = paper; circles = computed points; uncomputed radii remain empty');fig.savefig(OUT/'figure10_comparison.png',dpi=170);plt.close(fig)
 save(OUT/'comparison.json',{'fig9':comp9,'fig10':comp10,'normalization':'effective orbital energy flux per epsilon^2 zeta^2; unit mass source -> flux / alpha^6. No r0 ut adjustment.','reference_sha256':sha(OUT/'reference.json'),'limitations':['Finite numerical settings, no full convergence certification.','Sparse computed orbit samples are not a continuous solved curve.','Reference interpolation near zero/threshold is not a reliable relative-error bound.']})

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
  metric,audit=precompute_metric(DenseLorenzMetricMode(r,.88,mg,20),radii,theta,ROOT/'outputs/metric_cache',workers=workers)
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
  subprocess.run([sys.executable,str(OUT/'plot_factor_control.py')],stdout=subprocess.DEVNULL,check=True)
  subprocess.run([sys.executable,str(OUT/'plot_leaver_control.py')],stdout=subprocess.DEVNULL,check=True)

def run_all(workers):
 manifest=OUT/'execution.json';fp=fingerprint();plan=[(20.,mg) for mg in range(1,5)]+[(r,mg) for r in ORBIT_PLAN[1:] for mg in range(1,5)]+[(20.,mg) for mg in range(5,10)]
 if manifest.exists():
  previous=json.loads(manifest.read_text())
  if previous['implementation_sha256']!=fp:raise RuntimeError('Refuse resume with changed numerical implementation')
 for file in OUT.glob('r*_c*_l*_m*.json'):check_hashes(json.loads(file.read_text())['implementation_sha256'])
 state={'status':'running','pid':os.getpid(),'started_utc':stamp(),'groups':plan,'completed_groups':[],'failures':[],'implementation_sha256':fp,'orbit_plan':ORBIT_PLAN,'workers':workers,'limitations':['Finite truncation L20 -> j18, qmetric40, qsource64, pmax12.','Fig10 is a sparse orbit scan with extra samples near 18.3 and the dipole threshold.','Weak high-ell fluxes need precision/convergence checks.','Cloud BL normalization and complex-frequency policy remain explicit local conventions.']};save(manifest,state)
 for r,mg in plan:
  if fp!=fingerprint():raise RuntimeError('Numerical implementation changed; stop before mixing versions')
  state.update(active_orbit=r,active_metric_abs_m=mg,updated_utc=stamp());save(manifest,state)
  with (OUT/f'group_r{r:g}_mg{mg}.log').open('a') as log:
   proc=subprocess.run([sys.executable,'-u',__file__,'--group',str(r),str(mg),'--workers',str(workers)],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
  if proc.returncode:state['failures'].append({'r0':r,'mg':mg,'exit_code':proc.returncode})
  else:state['completed_groups'].append([r,mg])
  state['updated_utc']=stamp();save(manifest,state);render()
 state['status']='completed_finite_scan_convergence_pending' if not state['failures'] else 'partial_with_failures';state['completed_utc']=stamp();save(manifest,state)

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--run',action='store_true');ap.add_argument('--group',nargs=2,type=float);ap.add_argument('--workers',type=int,default=3);args=ap.parse_args()
 if args.group:group(args.group[0],int(args.group[1]),args.workers)
 else:
  reference();ingest();render()
  if args.run:
   try:run_all(args.workers)
   except Exception:
    p=OUT/'execution.json';d=json.loads(p.read_text());d.update(status='failed',traceback=traceback.format_exc());save(p,d);raise
