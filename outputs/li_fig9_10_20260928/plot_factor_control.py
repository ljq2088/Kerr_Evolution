"""Predetermined r0*u^t squared diagnostic; never changes a physical result."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=Path(__file__).resolve().parent
r=json.loads((OUT/'reference.json').read_text());d=json.loads((OUT/'comparison.json').read_text())
def factor(x):
 x=np.asarray(x);ut=(x**1.5+.88)/np.sqrt(x**3-3*x*x+2*.88*x**1.5);return (x*ut)**2
rows=[]
for row in d['fig10']:
 v=abs(row['local_signed']);f=float(factor(row['r0']));rows.append(dict(row,factor_r0ut_squared=f,local_times_factor=v*f,ratio_times_factor=v*f/row['Li_magnitude'],paper_over_local=row['Li_magnitude']/v if v else None))
cross=[]
for ref in r['fig10']:
 if ref['cloud']!=1 or ref['boundary']!='infinity':continue
 match=next(x for x in r['fig9'] if (x['ell'],x['m'])==(ref['ell'],ref['m']))
 value=float(10**np.interp(20.,ref['r'],np.log10(ref['Li_magnitude'])))
 cross.append(dict(ell=ref['ell'],m=ref['m'],Fig9=match['Li'],Fig10=value,Fig10_over_Fig9=value/match['Li'],ratio_over_r0ut_squared=value/match['Li']/float(factor(20.))))
fig,axs=plt.subplots(2,2,figsize=(13,8),layout='constrained',sharex='col',gridspec_kw={'height_ratios':[3,1]})
for c in (1,2):
 col=c-1
 for ref in r['fig10']:
  if ref['cloud']!=c:continue
  b,l,m=ref['boundary'],ref['ell'],ref['m'];style=['-','--','-.',':'][[v for v in r['fig10'] if v['cloud']==c and v['boundary']==b].index(ref)]
  color='C0' if b=='infinity' else 'C1';axs[0,col].semilogy(ref['r'],ref['Li_magnitude'],style,color=color,label=f'Li {b} ({l},{m})',lw=1.1)
  values=[v for v in rows if (v['cloud'],v['boundary'],v['ell'],v['m'])==(c,b,l,m)]
  for v in values:
   if v['local_times_factor']>0:axs[0,col].semilogy(v['r0'],v['local_times_factor'],'o',mfc='white',ms=4,color=color)
   if v['ratio_times_factor']>0:axs[1,col].semilogy(v['r0'],v['ratio_times_factor'],'o',color=color,ms=4)
 axs[0,col].set(title=f'Cloud ell_c=m_c={c}',ylabel='|Effective energy flux|',xlim=(4.1,50.1));axs[0,col].legend(fontsize=7,ncol=2);axs[0,col].grid(alpha=.2)
 axs[1,col].axhline(1,color='k',ls='--');axs[1,col].set(xlabel='r0/M',ylabel='Scaled local / Li');axs[1,col].grid(alpha=.2)
fig.suptitle('Hypothesis only: local flux multiplied by (r0 ut)^2\nLines: paper; circles: available computed points. Physical flux remains unchanged.')
fig.savefig(OUT/'figure10_factor_hypothesis.png',dpi=170);plt.close(fig)
(OUT/'factor_diagnostic.json').write_text(json.dumps({'rows':rows,'published_cross_figure_comparison':cross,'limitation':'A numerical scale match does not establish an unpublished author normalization or justify modifying physical flux.'},indent=2))
print(json.dumps({'published_cross_figure_comparison':cross,'cloud2_rows':[{k:v for k,v in row.items() if k in ('ell','m','boundary','paper_over_local','ratio_times_factor')} for row in rows if row['cloud']==2]},indent=2))
