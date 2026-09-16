"""Compare unchanged finite fluxes with Li et al's published comparison curves.

Uses original vector PDFs, not reconstructed author raw numerical arrays.
No normalization is fitted or applied to local solver output.
"""
import hashlib,json
from pathlib import Path
import numpy as np
import pymupdf as fitz
from scipy.interpolate import PchipInterpolator


def points(d):
    assert all(i[0]=='l' for i in d['items'])
    return np.array([tuple(d['items'][0][1])]+[tuple(i[2]) for i in d['items']])


def main():
    root=Path(__file__).resolve().parents[1];out=root/'docs/environment_reproduction'
    oldpdf=root/'outputs/environment_reference/Flux_Inf_Hor_TwoPanels_2.pdf'
    newpdf=root/'outputs/paper_original_reference/li_2507_02045v2/source/total_flux_11.pdf'
    old=fitz.open(oldpdf)[0].get_drawings();new=fitz.open(newpdf)[0].get_drawings()
    grids=[d for d in old if d['color'] and max(abs(c-.75) for c in d['color'])<1e-3 and d['rect'].x0>230]
    xs=sorted(a.x for d in grids for kind,a,b in d['items'] if kind=='l' and a.x==b.x)
    ys=sorted(a.y for d in grids for kind,a,b in d['items'] if kind=='l' and a.y==b.y)
    assert len(xs)==5 and len(ys)==4
    oldx=np.polyfit(xs,[10,20,30,40,50],1);oldy=np.polyfit(ys,[-1,-2,-3,-4],1)
    oldcurves={}
    for d in old:
        if d['color'] and max(abs(c-v) for c,v in zip(d['color'],[.88072,.61104,.14205]))<1e-5 and d['rect'].x0>230 and len(d['items'])>=20:
            pts=points(d);oldcurves['infinity' if d['dashes']=='[] 0' else 'horizon']=np.column_stack((np.polyval(oldx,pts[:,0]),np.polyval(oldy,pts[:,1])))
    nx=np.polyfit([61.12561798095703,98.98464965820312,136.8436737060547,174.7027130126953,212.56173706054688],[10,20,30,40,50],1)
    ny=np.polyfit([148.2509307861328,110.11843872070312,71.98595428466797,33.85345458984375],[-4,-3,-2,-1],1)
    curves={}
    for i,name in [(41,'li_horizon'),(42,'li_infinity'),(43,'dyson_new_horizon'),(44,'dyson_new_infinity')]:
        pts=points(new[i]);curves[name]=np.column_stack((np.polyval(nx,pts[:,0]),np.polyval(ny,pts[:,1])))
    def at(curve,r):return float(10**np.interp(r,curve[:,0],curve[:,1])*.3**6)
    def interp_indicator(curve,r):return float(10**PchipInterpolator(curve[:,0],curve[:,1])(r)*.3**6/at(curve,r)-1)
    localpath=out/'figure2_three_orbit_comparison.json';local=json.loads(localpath.read_text());rows=[]
    factor=2/(1+np.sqrt(1-.88**2))
    for row in local['rows']:
        rp=row['rp']
        for boundary in ['horizon','infinity']:
            original=at(oldcurves[boundary],rp);updated=at(curves['dyson_new_'+boundary],rp);li=at(curves['li_'+boundary],rp);computed=abs(row['fluxes'][boundary]['computed_signed'])
            rows.append(dict(rp=rp,boundary=boundary,local_finite_magnitude=computed,dyson_original_magnitude=original,dyson_in_later_comparison_magnitude=updated,li_magnitude=li,new_over_original=updated/original,local_vs_original_percent=100*(computed/original-1),local_vs_updated_dyson_percent=100*(computed/updated-1),local_vs_li_percent=100*(computed/li-1),new_dyson_log_pchip_vs_linear_relative=interp_indicator(curves['dyson_new_'+boundary],rp),li_log_pchip_vs_linear_relative=interp_indicator(curves['li_'+boundary],rp)))
    coverage_hashes={}
    for row in rows:
        suffix='' if row['rp']==20 else '_rp'+str(row['rp'])
        cp=out/f'flux_coverage_L18_nt18_i6_h5_f5{suffix}.json'
        coverage=json.loads(cp.read_text());coverage_hashes[str(cp.relative_to(root))]=hashlib.sha256(cp.read_bytes()).hexdigest()
        block=coverage[row['boundary']]
        total=sum(c['flux'][row['boundary']]['orbital_energy'] for c in block['modes'] if c['ell']<=5)
        row['local_same_ellmax5_magnitude']=abs(total)
        row['local_same_ellmax5_vs_li_percent']=100*(abs(total)/row['li_magnitude']-1)
        row['local_same_ellmax5_channels']=sum(c['ell']<=5 for c in block['modes'])
    scan=[]
    for r in np.linspace(5,40,141):
        rh=at(curves['dyson_new_horizon'],r)/at(oldcurves['horizon'],r);ri=at(curves['dyson_new_infinity'],r)/at(oldcurves['infinity'],r)
        scan.append(dict(rp=float(r),horizon_new_over_old=rh,infinity_new_over_old=ri,horizon_ratio_over_area_factor=rh/factor))
    result=dict(status='published_reference_normalization_change_explains_bulk_offset_residual_not_closed',references=dict(old='https://arxiv.org/src/2501.09806v1',new='https://arxiv.org/src/2507.02045v2',new_footnote='Footnote1 explicitly identifies scalar separation-constant and horizon-flux normalization errors in Dyson et al.; no old-code expression is given.'),inputs={str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [oldpdf,newpdf,localpath,Path(__file__)]},same_cutoff_coverage_sha256=coverage_hashes,normalization='per q^2 eta; eta=M_cloud/M. Plot normalization q^2 epsilon^2 converted by alpha^6, alpha=.3.',rows=rows,area_factor_hypothesis=dict(formula='(r_plus^2+a^2)/r_plus^2=2M/r_plus',value_at_label_spin_088=factor,value_at_exact_threshold_spin=2/(1+np.sqrt(1-.8771530275949366**2)),status='Geometric factor and quantitative PDF fingerprint, NOT the unpublished original code expression.',scan_radius_range=[5,40],max_horizon_fractional_departure=float(max(abs(v['horizon_ratio_over_area_factor']-1) for v in scan)),max_infinity_new_old_departure=float(max(abs(v['infinity_new_over_old']-1) for v in scan))),scan=scan,limitations=['PDF polylines are digitized published curves, not author raw arrays; no rigorous readout error bound.','Li actually uses a=.88, while local/Dyson exact threshold spin is about.877153; comparisons are not identical parameter runs.','Local totals retain historical finite source and mode cutoffs.','Separation-constant correction is explicitly mentioned in the paper, but its exact old implementation is not supplied.','Physical production source, radial solver and flux functions are unchanged.','Li total curves use scalar ell<=5 at both boundaries; same_ellmax5 comparisons explicitly match this cutoff. Existing local infinity totals use ell<=6.'])
    path=out/'li_reference_comparison_20260916.json';path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10})
    fig,axes=plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
    rr=np.linspace(8,32,400)
    for c,label,col,style in [(oldcurves['horizon'],'Dyson original Fig. 2','#777777','--'),(curves['dyson_new_horizon'],'Dyson curve in Li Fig. 2','#D55E00','-'),(curves['li_horizon'],'Li et al. (a/M = 0.88)','#009E73','-')]:
        axes[0].semilogy(rr,[at(c,r) for r in rr],style,color=col,label=label)
    hs=[r for r in rows if r['boundary']=='horizon'];rx=[r['rp'] for r in hs];fy=[r['local_finite_magnitude'] for r in hs]
    axes[0].semilogy(rx,fy,'o',color='#0072B2',markersize=6,label='Local unchanged finite sums',zorder=5)
    axes[0].set(xlabel=r'$r_p/M$',ylabel=r'$|F^{s,H}|/[q^2(M_c/M)]$',title='Horizon flux: original and later references',xlim=(8,32),xticks=[10,20,30]);axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
    for key,label,col,style in [('local_vs_original_percent','vs. original Dyson','#777777','--'),('local_vs_updated_dyson_percent','vs. Dyson in Li','#D55E00','-'),('local_vs_li_percent','vs. Li (different spin)','#009E73','-')]:
        axes[1].plot(rx,[r[key] for r in hs],style,marker='o',color=col,label=label)
    axes[1].axhline(0,color='black',lw=.7);axes[1].set(xlabel=r'$r_p/M$',ylabel='Local / reference - 1  [%]',title='Large offset becomes a few-percent residual',xticks=[10,20,30],xlim=(8,32));axes[1].grid(alpha=.2);axes[1].legend(fontsize=8)
    fig.suptitle(r'$\alpha=0.3$: published reference change; no fitted rescaling of local flux',fontsize=12)
    fig.savefig(out/'li_reference_comparison_20260916.png',dpi=180);plt.close(fig)
    print(json.dumps(dict(rows=rows,area_factor_hypothesis=result['area_factor_hypothesis']),indent=2))


if __name__=='__main__':main()
