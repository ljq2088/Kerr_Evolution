"""Compare cached |Sigma-source projections| with Li/Dyson vector curves.

No source solve, amplitude fitting, or normalization inferred from numbers.
The literal alpha^-3 follows the stated amplitude expansion convention.
Its identification with the plotted source units is UNRESOLVED, not assumed.
"""
from pathlib import Path
import csv,hashlib,json
import numpy as np
from numpy.polynomial.legendre import legfit,legval
from scipy.interpolate import PchipInterpolator
import pymupdf as fitz
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"docs/environment_reproduction"
PDF=ROOT/"outputs/paper_original_reference/li_2507_02045v2/source/source_l2l3.pdf"
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
NAMES={2:'forced_mode_nr8_nt18_L18_mg1_sl2_gh0.0001_go4000_inner0.0005_outer320_log_h32.json',3:'forced_mode_nr8_nt18_L18_inner0.0005_outer320_log_h32.json'}

def cloud_normalization_check():
 # Independent raw Leaver a0=1 is a explicitly specified mathematical check;
 # it is not asserted to be the unpublished source-plot convention.
 from environment_source import ThresholdCloud
 from paper_leaver_cloud import LeaverThresholdCloud
 c=ThresholdCloud();l=LeaverThresholdCloud(terms=200,dps=40)
 mass,charge=l.integrals();raw_h=float(np.exp(-float(l.k*l.rp))*float(l.rp-l.rm)**float(l.beta))
 cp=ROOT/'outputs/paper_original_reference/li_2507_02045v2/source/scalar_profile.pdf'
 page=fitz.open(cp)[0];items=page.get_drawings()[2]['items'];points=np.array([tuple(items[0][1])]+[tuple(z[2]) for z in items])
 x=(points[:,0]-36.552669525146484)*250/(228.3687286376953-36.552669525146484)
 y=(140.0643768310547-points[:,1])*.4/(140.0643768310547-20.223236083984375)
 samples=[]
 for r in [5.,10.,20.,50.,100.]:
  local=float(c.radial(r)[0,0]/c.mu**3);paper=float(np.interp(r,x,y))
  samples.append(dict(r=r,local_R_literal_zeta_units=local,paper_R_pdf_polyline=paper,local_over_paper=local/paper,raw_a0_1_Leaver_R=float(l.radial(r)[0])))
 return dict(status='cloud_profile_is_consistent_with_literal_zeta_units_but_source_figure_is_not',
  cloud_profile_pdf_sha256=sha(cp),raw_Leaver_ansatz='exp(-k*r)*(r-r_minus)^beta*sum a_n*((r-r_plus)/(r-r_minus))^n with a0=1',
  raw_Leaver_unit_a0_cloud_mass=float(mass),raw_Leaver_horizon_R=raw_h,
  raw_Leaver_to_literal_zeta_amplitude=float(1/np.sqrt(mass)/c.mu**3),
  horizon_unit_to_unit_mass_amplitude=float(c.amplitude),horizon_unit_to_literal_zeta_amplitude=float(c.amplitude/c.mu**3),
  samples=samples,limitation='Cloud a is synchronized locally vs 0.88 in reference. Raw a0=1 convention is not identified as the author source plot convention.')

def main():
 doc=fitz.open(PDF);draw=doc[0].get_drawings()
 configs={2:dict(curves={'Li':2,'Dyson_replotted_by_Li':3},xticks=[105.03203582763672,174.45848083496094],yticks=[139.0921173095703,120.44125366210938,101.79039764404297,83.13954162597656,64.48867797851562,45.83781433105469,27.186965942382812],ylog=list(range(-7,0))),
 3:dict(curves={'Li':29,'Dyson_replotted_by_Li':30},xticks=[352.7120361328125,422.13848876953125],yticks=[135.60842895507812,115.2350082397461,94.8615951538086,74.4881820678711,54.114768981933594,33.741363525390625],ylog=list(range(-7,-1)))}
 curves={};allvertices=[];rows=[];node_rows=[];inputs=[PDF,Path(__file__)];params={}
 for ell,cfg in configs.items():
  cfg['xmap']=np.polyfit(cfg['xticks'],[1,2],1).tolist();cfg['ymap']=np.polyfit(cfg['yticks'],cfg['ylog'],1).tolist()
  for name,idx in cfg['curves'].items():
   items=draw[idx]['items'];assert all(z[0]=='l' for z in items)
   pts=np.array([tuple(items[0][1])]+[tuple(z[2]) for z in items]);assert np.all(np.diff(pts[:,0])>0)
   logr=np.polyval(cfg['xmap'],pts[:,0]);logv=np.polyval(cfg['ymap'],pts[:,1]);curves[(ell,name)]=(logr,logv)
   allvertices.extend(dict(ell=ell,m=ell,reference=name,r=float(10**rr),abs_source=float(10**ss),pdf_x=float(x),pdf_y=float(y)) for (rr,ss,(x,y)) in zip(logr,logv,pts))
  p=OUT/NAMES[ell];inputs.append(p);data=json.loads(p.read_text());par=data['parameters'];params[str(ell)]=par
  assert par['scalar_ell']==par['scalar_m']==ell and par['cloud_mass']==1 and par['metric']['orbital_radius']==20
  rad=np.array([x['r'] for x in data['samples']]);source=np.array([complex(*x['source']) for x in data['samples']]);panels=par['source_panels'];a=par['metric']['a'];rplus=1+np.sqrt(1-a*a)
  scale=par['alpha']**-3/np.sqrt(par['cloud_mass']);parts=[]
  for i,(lo,hi) in enumerate(zip(panels[:-1],panels[1:])):
   mask=(rad>lo)&(rad<hi);rr=rad[mask];ss=source[mask];islog=(i==0 and par.get('horizon_log_first_panel',False));z=np.log(rr-rplus) if islog else rr;zlo,zhi=np.log([lo-rplus,hi-rplus]) if islog else (lo,hi)
   x=2*(z-zlo)/(zhi-zlo)-1;assert np.allclose(x,np.polynomial.legendre.leggauss(len(x))[0],rtol=0,atol=2e-9)
   parts.append(dict(lo=lo,hi=hi,zlo=zlo,zhi=zhi,log=islog,coeff=legfit(x,ss,len(x)-1),real_pchip=PchipInterpolator(z,ss.real),imag_pchip=PchipInterpolator(z,ss.imag),radii=rr))
  def evaluate(r):
   part=next(p for p in parts if p['lo']<r<p['hi']);z=np.log(r-rplus) if part['log'] else r;x=2*(z-part['zlo'])/(part['zhi']-part['zlo'])-1
   val=legval(x,part['coeff']);alt=part['real_pchip'](z)+1j*part['imag_pchip'](z)
   return val,alt,part
  for target in [3.5,5.,8.,10.,15.,25.,50.,100.,200.]:
   value,alt,part=evaluate(target);fixed=dict(ell=ell,m=ell,r=target,local_source_real=float(value.real),local_source_imag=float(value.imag),
     cloud_amplitude_conversion=scale,local_abs_source_literal_zeta_units=float(abs(value)*scale),
     local_interpolation='same panel degree N-1 complex Legendre; log r-r+ on first panel',
     panel_bounds=[part['lo'],part['hi']],interpolation_method_change=float(abs(alt-value)/(abs(value) or 1)))
   j=int(np.argmin(abs(rad-target)));raw=dict(ell=ell,m=ell,requested_near_r=target,r=float(rad[j]),
      local_abs_source_literal_zeta_units=float(abs(source[j])*scale),local_interpolation='none: original cached Gauss sample')
   for row in [fixed,raw]:
    for name in cfg['curves']:
     x,y=curves[(ell,name)];reference=float(10**np.interp(np.log10(row['r']),x,y));row[name+'_abs_source']=reference;row['local_over_'+name]=row['local_abs_source_literal_zeta_units']/reference
    row['Li_over_Dyson']=row['Li_abs_source']/row['Dyson_replotted_by_Li_abs_source']
   rows.append(fixed);node_rows.append(raw)
 # This dimensionless shape ratio does not alter either absolute data set.
 for ell in [2,3]:
  selected=[x for x in rows if x['ell']==ell]
  anchor=next(x for x in selected if x['r']==10.)
  for row in selected:
   for reference in ['Li','Dyson_replotted_by_Li']:
    row['shape_ratio_relative_to_r10_'+reference]=row['local_over_'+reference]/anchor['local_over_'+reference]
 cloud_info=cloud_normalization_check()
 with (OUT/'li_source22_33_vector_vertices_20260917.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(allvertices[0]));w.writeheader();w.writerows(allvertices)
 for name,records in [('li_source22_33_comparison_20260917',rows),('li_source22_33_raw_nodes_20260917',node_rows)]:
  with (OUT/f'{name}.csv').open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
 fig,axes=plt.subplots(2,2,figsize=(11.4,7.1),sharex='col',gridspec_kw={'height_ratios':[2,1]},layout='constrained')
 for col,ell in enumerate([2,3]):
  for name,color,style in [('Li','#1764ab','-'),('Dyson_replotted_by_Li','#c44202','--')]:
   x,y=curves[(ell,name)];axes[0,col].loglog(10**x,10**y,color=color,ls=style,lw=1.7,label=name.replace('_',' '))
  rr=[x for x in rows if x['ell']==ell];xx=[x['r'] for x in rr]
  axes[0,col].loglog(xx,[x['local_abs_source_literal_zeta_units'] for x in rr],'o',color='#222222',ms=4,label='Local: fixed radii')
  for name,color in [('Li','#1764ab'),('Dyson_replotted_by_Li','#c44202')]:axes[1,col].semilogx(xx,[x['local_over_'+name] for x in rr],'o-',color=color,ms=3,label=name)
  axes[0,col].set(title=fr'Scalar $(\ell,m)=({ell},{ell})$, $r_p=20M$',ylabel=r'$|\mathscr{S}_{\ell m}|$')
  axes[1,col].axhline(1,color='#555555',lw=.8);axes[1,col].set(xlim=(3,250),xlabel='r/M',ylabel='Local / reference')
  axes[0,col].legend(fontsize=8);axes[1,col].grid(alpha=.2)
 fig.suptitle('Source magnitudes: plotted absolute normalization remains unresolved',fontsize=13)
 fig.supxlabel(r'Local $J$ multiplied only by $\alpha^{-3}$ under the literal expansion convention; no fitted correction. Spins differ: 0.877153 vs 0.88.',fontsize=9)
 figure=OUT/'li_source22_33_comparison_20260917.png';fig.savefig(figure,dpi=180);plt.close(fig)
 report=dict(status='two_mode_radial_source_comparison_vector_digitization_not_author_arrays',source_url='https://arxiv.org/html/2507.02045v2',paper_source_pdf=str(PDF.relative_to(ROOT)),
  exact_variable_relation=dict(local='J_lm=integral S_lm* Sigma h^{ab} Hessian_ab(phi_cloud) dOmega',
   Li='mathscr S_lm in d_r(Delta d_r R)+V R=mathscr S',
   equation_relation='For the same normalized h and physical cloud, the coefficient in the literal expansion Phi=zeta Phi10+epsilon*zeta Phi11 is J_local/(alpha^3 sqrt(eta)); identification with the plotted source normalization is unresolved',
   extra_Delta_factor=False,extra_sqrt_r2_plus_a2_factor=False,angular_normalization='integral |S(theta)|^2 sin(theta) dtheta=1/(2*pi)',
   tex_lines=[148,164,371,373,438,448,491,511,851]),
  paper_parameters=dict(a=.88,alpha=.3,rp=20,cloud_ell=1,cloud_m=1,cloud_n=0),local_parameters=params,
  normalization_status='UNRESOLVED: literal zeta-stripped source is about 21.6 times the plotted source; no fitted factor was applied',
  cloud_normalization_check=cloud_info,curve_calibration=configs,rows=rows,raw_node_rows=node_rows,
  inputs_sha256={str(p.relative_to(ROOT)):sha(p) for p in inputs},
  outputs_sha256={str(p.relative_to(ROOT)):sha(p) for p in [figure,OUT/'li_source22_33_vector_vertices_20260917.csv',OUT/'li_source22_33_comparison_20260917.csv',OUT/'li_source22_33_raw_nodes_20260917.csv']},
  limitations=['Absolute source-plot normalization is not established by the public plotting code, which is unavailable; an amplitude mismatch is not proof of a source error.',
   'Ratios use the literal alpha^-3 amplitude convention. Shape ratios relative to r10 are separately labeled and never applied to the absolute data.',
   'Reference is the absolute value of a vector curve, not complex source samples; phase is untested.',
   'Reference interpolation is exactly linear on the log-log PDF polyline, not assumed true author data between vertices.',
   'Local polynomial interpolation follows recorded source panels; PCHIP change is an interpolation diagnostic, not a rigorous solver error.',
   'Original Gauss-node comparison avoids local interpolation and is included separately.',
   'Spin/cloud frequencies differ between local synchronized a=.877153... and Li a=.88; Dyson source is the curve replotted in Li, not an independently downloaded author array.',
   'Metric/angular/source convergence is not inferred from agreement of these two low-order source magnitudes.'])
 (OUT/'li_source22_33_comparison_20260917.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
 print(json.dumps({'rows':rows,'nodes':node_rows},indent=2))

if __name__=='__main__':main()
