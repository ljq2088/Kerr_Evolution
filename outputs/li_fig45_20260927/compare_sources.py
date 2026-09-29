"""Li Fig4/5 source comparison; raw physical sources retained without rescaling."""
import json,sys,hashlib,os,traceback
from pathlib import Path
import numpy as np
import pymupdf as fitz
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
BANK=ROOT/'docs/li_alignment/rp20_cloud11_j18_L20_q40_nr12_h64_cache512'
sys.path.insert(0,str(ROOT/'src'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
 t=p.with_suffix('.tmp');t.write_text(json.dumps(d,indent=2));t.replace(p)
def refs():
 result={}
 configs=[('source_l2l3',[(2,[2],[2], [105.03203582763672,174.45848083496094],[139.0921173095703,27.186965942382812],[-7,-1]),(3,[3],[29],[352.7120361328125,422.13848876953125],[135.60842895507812,33.741363525390625],[-7,-2])]),('source_l6l7',[(6,[6,4,2],[2,4,6],[105.75375366210938,172.35755920410156],[134.4059295654297,26.999130249023438],[-12,-2]),(7,[7,5,3],[31,33,35],[353.4337463378906,420.0375671386719],[131.37142944335938,26.385711669921875],[-12,-2])])]
 for name,panels in configs:
  path=ROOT/f'outputs/paper_original_reference/li_2507_02045v2/source/{name}.pdf';draw=fitz.open(path)[0].get_drawings()
  for l,ms,ids,xt,yt,yl in panels:
   for m,i in zip(ms,ids):
    row={'ell':l,'m':m,'pdf_sha256':sha(path)}
    for who,j in [('Li',i),('Dyson',i+1)]:
     items=draw[j]['items'];assert all(v[0]=='l' for v in items)
     v=np.array([tuple(items[0][1])]+[tuple(z[2]) for z in items]);assert np.all(np.diff(v[:,0])>0)
     row[who]={'r':(10**np.polyval(np.polyfit(xt,[1,2],1),v[:,0])).tolist(),'source':(10**np.polyval(np.polyfit(yt,yl,1),v[:,1])).tolist()}
    result[f'{l}_{m}']=row
 save(OUT/'reference_curves.json',result);return result

def render():
 reference=refs();a=.88;r0=20.;ut=(r0**1.5+a)/np.sqrt(r0**3-3*r0*r0+2*a*r0**1.5);f=r0*ut
 summary={'ut':ut,'r0_ut':f,'normalization':'physical unit-cloud-mass source / alpha^3; no fitted factor','modes':[],'limitations':['PDF vector digitization; no rigorous error bars.','r0*u^t curves are explicit hypotheses, not solver corrections.','Finite L20 -> spherical j18, angular q64, pmax12; convergence not certified.']}
 for fig,ells in [(4,[2,3]),(5,[6,7])]:
  image,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained',sharex='col')
  for col,l in enumerate(ells):
   for k,ref in reference.items():
    if ref['ell']!=l:continue
    m=ref['m'];color={2:'C0',3:'C0',4:'C1',5:'C1',6:'C2',7:'C2'}[m]
    rr=np.array(ref['Li']['r']);ss=np.array(ref['Li']['source'])
    axes[0,col].loglog(rr,ss,color=color,label=f'Li m={m}')
    axes[0,col].loglog(ref['Dyson']['r'],ref['Dyson']['source'],'--',color=color,alpha=.5)
    path=OUT/f'source_l{l}_m{m}.json'
    if not path.exists():continue
    d=json.loads(path.read_text());r=np.array(d['r']);z=np.array(d['source']);v=abs(z[:,0]+1j*z[:,1])/.3**3
    mask=(r>=rr.min())&(r<=rr.max());r=r[mask];v=v[mask]
    paper=10**np.interp(np.log10(r),np.log10(rr),np.log10(ss));ratio=v/paper
    axes[0,col].loglog(r,v,'.',color=color,label=f'Local raw m={m}')
    axes[1,col].semilogx(r,ratio,'o-',ms=2,color=color,label=f'raw m={m}')
    axes[1,col].semilogx(r,ratio/f,':',color=color,label=f'raw / (r0 ut), m={m}')
    central=(r>=3)&(r<=250)&(abs(r-20)>1)
    summary['modes'].append({'ell':l,'m':m,'origin':d['origin'],'input_sha256':sha(path),'raw_ratio_median':float(np.median(ratio[central])),'raw_ratio_minmax':[float(ratio[central].min()),float(ratio[central].max())],'divided_r0ut_ratio_median':float(np.median(ratio[central]/f)),'points':[{'r':float(x),'local_over_Li':float(y),'divided_r0ut_over_Li':float(y/f)} for x,y in zip(r,ratio)]})
   axes[0,col].set(title=f'ell={l}',ylabel='|source| (Li expansion units)');axes[0,col].legend(fontsize=8);axes[0,col].grid(alpha=.2)
   axes[1,col].axhline(1,color='k',ls='--');axes[1,col].set(xlabel='r/M',ylabel='Local / Li',yscale='log');axes[1,col].legend(fontsize=7);axes[1,col].grid(alpha=.2)
  image.suptitle(f'Li Fig.{fig}: absolute source comparison (dashed: Dyson replotted by Li)')
  image.savefig(OUT/f'figure{fig}_comparison.png',dpi=160);plt.close(image)
 save(OUT/'comparison.json',summary)

def ingest_existing():
 for l,m in [(2,2),(3,3)]:
  p=BANK/f'mode_l{l}_m{m}.json';d=json.loads(p.read_text())
  for module,digest in d['implementation_sha256'].items():
   src=ROOT/'src'/(module if module.endswith('.json') else module+'.py')
   if not src.exists() or sha(src)!=digest:raise ValueError(f'Implementation mismatch: {module}')
  assert d['config']['a']==.88 and d['config']['r0']==20. and d['config']['cloud_ell']==1
  save(OUT/f'source_l{l}_m{m}.json',{'origin':'reused completed Li-aligned calculation; all recorded code hashes match','input':str(p),'input_sha256':sha(p),'r':d['source_radii'],'source':d['source']})

def fresh():
 from li_normalized_cloud import LiNormalizedCloud
 from li_separated_source import LiSeparatedSource
 from environment_source import angular_mode
 from report_li_aligned_flux import spherical_coefficients
 from environment_metric_sampling import precompute_metric
 from environment_dense_metric import DenseLorenzMetricMode
 os.nice(10);cloud=LiNormalizedCloud()
 status={'status':'running','pid':os.getpid(),'completed':[],'expected':[[6,2],[7,3],[6,4],[7,5],[6,6],[7,7]]};save(OUT/'execution.json',status)
 full=np.array(json.loads((BANK/'mode_l2_m2.json').read_text())['source_radii']);radii=full[(full>=3)&(full<=300)]
 x,w=np.polynomial.legendre.leggauss(40);theta=np.arccos(x)
 for l,m in status['expected']:
  mg=m-1;path=OUT/f'source_l{l}_m{m}.json'
  if path.exists():status['completed'].append([l,m]);continue
  bp=BANK/f'metric_spherical_mg{mg}.npz'
  if bp.exists():
   with np.load(bp,allow_pickle=False) as b:
    meta=json.loads(str(b['metadata']));allr=b['radii'];coef=b['coefficients']
   for module,digest in meta['implementation_sha256'].items():
    src=ROOT/'src'/(module if module.endswith('.json') else module+'.py')
    if sha(src)!=digest:raise ValueError('Saved bank source code changed')
   coefficients=np.array([coef[np.flatnonzero(allr==r)[0]] for r in radii]);audit={'bank':str(bp),'sha256':sha(bp)}
  else:
   metric,audit=precompute_metric(DenseLorenzMetricMode(20.,.88,mg,20),radii,theta,ROOT/'outputs/metric_cache',workers=1)
   coefficients=np.array([spherical_coefficients(metric,r,theta,w,.88,mg,18) for r in radii])
  omega=cloud.omega+mg/(20**1.5+.88)
  projector=LiSeparatedSource(a=.88,omega_c=cloud.omega,m_c=1,metric_m=mg,spherical_lmax=18,cloud_angular=cloud.angular_state,target_angular=lambda t:angular_mode(t,l,m,.88**2*(omega**2-.3**2))[0],quadrature=64,pmax=12,dps=64)
  result=[]
  for index,(r,c) in enumerate(zip(radii,coefficients)):
   z=projector.project(r,cloud.radial_state(r)[:,0],c);result.append([float(z.real),float(z.imag)])
   status.update(active_mode=[l,m],source_points_done=index+1,source_points_total=len(radii));save(OUT/'execution.json',status)
   save(OUT/f'checkpoint_l{l}_m{m}.json',{'r':radii[:index+1].tolist(),'source':result,'metric':audit})
  save(path,{'origin':'new separated-source projection from provenance-checked metric data','metric':audit,'r':radii.tolist(),'source':result,'cloud':cloud.provenance})
  status['completed'].append([l,m]);save(OUT/'execution.json',status);render()
 status['status']='completed_finite_resolution_comparison';save(OUT/'execution.json',status)
if __name__=='__main__':
 try:
  ingest_existing();render()
  if '--run' in sys.argv:fresh()
 except Exception:
  save(OUT/'failure.json',{'error':traceback.format_exc()});raise
