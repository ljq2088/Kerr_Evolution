"""Independent unshifted a0=1 Leaver normalization diagnostic; no fitting."""
from pathlib import Path
import sys,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1];sys.path.insert(0,str(ROOT/'src'))
from li_normalized_cloud import LiNormalizedCloud
factors={}
for c in [1,2]:
 cloud=LiNormalizedCloud(ell=c);conversion=np.exp(cloud.q*cloud.rp)*(cloud.rp-cloud.rm)**cloud.beta
 factors[c]={'raw_flux_scale':float(cloud.provenance['raw_bl_mass']*.3**6*abs(conversion)**2),'prefactor_conversion_abs':float(abs(conversion)),'cloud_provenance':cloud.provenance}
reference=json.loads((OUT/'reference.json').read_text());records=[json.loads(p.read_text()) for p in OUT.glob('r*_c*_l*_m*.json')];rows=[]
fig,axs=plt.subplots(2,2,figsize=(13,8),layout='constrained',sharex='col',gridspec_kw={'height_ratios':[3,1]})
for c in [1,2]:
 col=c-1;scale=factors[c]['raw_flux_scale']
 for ref in reference['fig10']:
  if ref['cloud']!=c:continue
  l,m,b=ref['ell'],ref['m'],ref['boundary'];j=[v for v in reference['fig10'] if v['cloud']==c and v['boundary']==b].index(ref);color='C0' if b=='infinity' else 'C1'
  axs[0,col].semilogy(ref['r'],ref['Li_magnitude'],['-','--','-.',':'][j],color=color,label=f'Li {b} ({l},{m})',lw=1.1)
  for d in records:
   if (d['cloud'],d['ell'],d['m'])!=(c,l,m) or not min(ref['r'])<=d['r0']<=max(ref['r']):continue
   value=abs(d['flux'][b])*scale;paper=float(10**np.interp(d['r0'],ref['r'],np.log10(ref['Li_magnitude'])));ratio=value/paper
   rows.append(dict(cloud=c,r0=d['r0'],ell=l,m=m,boundary=b,physical_flux=d['flux'][b],raw_Leaver_flux=value,Li_magnitude=paper,ratio=ratio))
   if value>0:axs[0,col].semilogy(d['r0'],value,'o',mfc='white',color=color,ms=4);axs[1,col].semilogy(d['r0'],ratio,'o',color=color,ms=4)
 axs[0,col].set(title=f'Cloud {c}: predicted scale = {scale:.7g}',ylabel='|Energy flux|',xlim=(4.1,50.1));axs[0,col].legend(fontsize=7,ncol=2);axs[0,col].grid(alpha=.2)
 axs[1,col].axhline(1,color='k',ls='--');axs[1,col].set(xlabel='r0/M',ylabel='Raw-Leaver local / Li');axs[1,col].grid(alpha=.2)
fig.suptitle('Fig.10 normalization hypothesis: unshifted Leaver radial series with a0=1\nLines: paper; circles: computed points. Independent constants; no fitted amplitude.')
fig.savefig(OUT/'figure10_raw_leaver_control.png',dpi=170);plt.close(fig)
result={'status':'normalization_hypothesis_not_author_code_confirmation','definition':'R_raw = exp(q*r)*(r-r_minus)^beta*x^s*sum a_n*x^n, a0=1; x=(r-r_plus)/(r-r_minus)','flux_scale_formula':'alpha^6 * raw_shifted_BL_mass * abs(exp(q*r_plus)*(r_plus-r_minus)^beta)^2','factors':factors,'rows':rows,'limitations':['No changes to physical sources or fluxes.','Constants depend on cloud state, not orbit radius. Second-orbit calculations distinguish this from (r0*u^t)^2.','Horizon channels have additional disagreements; no full reproduction claimed.']}
(OUT/'raw_leaver_normalization_control.json').write_text(json.dumps(result,indent=2))
print(json.dumps([r for r in rows if r['boundary']=='infinity'],indent=2))
