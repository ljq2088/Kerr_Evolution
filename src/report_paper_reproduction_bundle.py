"""Build a traceable scientific figure set from finite archived solver results.

No reference is used to renormalize a local solution. New runs are attached
separately because partial channels must never be drawn as total fluxes.
"""
import hashlib
import json
import shutil
from pathlib import Path
import numpy as np
import pymupdf as fitz
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Circle
from report_li_reference_comparison import points

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'docs/environment_reproduction'
OUT=ROOT/'docs/paper_reproduction_20260916'
INPUTS={}

def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name):
    p=DATA/name; INPUTS[str(p.relative_to(ROOT))]=digest(p)
    return json.loads(p.read_text())
def array(name):
    p=DATA/name; INPUTS[str(p.relative_to(ROOT))]=digest(p)
    return np.load(p)
def save(name,obj):
    (OUT/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def plot_save(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=190)
    fig.savefig(OUT/(name+'.svg'))
    svg=OUT/(name+'.svg'); svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    plt.close(fig)
def validate_modes(coverage):
    for boundary in ('infinity','horizon'):
        b=coverage[boundary]
        assert b['computed_count']==b['required_count']
        seen=set();values=[]
        for row in b['modes']:
            mode=(row['ell'],row['m']); assert mode not in seen;seen.add(mode)
            p=DATA/row['file']; actual=json.loads(p.read_text())
            assert actual['flux']==row['flux']
            assert (actual['parameters']['scalar_ell'],actual['parameters']['scalar_m'])==mode
            INPUTS[str(p.relative_to(ROOT))]=digest(p)
            values.append(row['flux'][boundary]['orbital_energy'])
        np.testing.assert_allclose(sum(values),b['finite_resolution_total'],rtol=1e-14,atol=0)

def flux_figure():
    comparison=load('li_reference_comparison_20260916.json')
    pdf=ROOT/'outputs/paper_original_reference/li_2507_02045v2/source/total_flux_11.pdf'
    INPUTS[str(pdf.relative_to(ROOT))]=digest(pdf)
    assert comparison['inputs'][str(pdf.relative_to(ROOT))]==digest(pdf)
    d=fitz.open(pdf); drawings=d[0].get_drawings()
    nx=np.polyfit([61.12561798095703,98.98464965820312,136.8436737060547,174.7027130126953,212.56173706054688],[10,20,30,40,50],1)
    ny=np.polyfit([148.2509307861328,110.11843872070312,71.98595428466797,33.85345458984375],[-4,-3,-2,-1],1)
    curves={}
    for i,key in [(41,'li_horizon'),(42,'li_infinity'),(43,'dyson_horizon'),(44,'dyson_infinity')]:
        p=points(drawings[i]); curves[key]=np.column_stack([np.polyval(nx,p[:,0]),10**np.polyval(ny,p[:,1])])
    rows=[]
    for rp,name in [(10,'flux_coverage_L18_nt18_i6_h5_f5_rp10.json'),(20,'flux_coverage_L18_nt18_i6_h5_f5_m0go4000_m2go4000.json'),(30,'flux_coverage_L18_nt18_i6_h5_f5_rp30.json')]:
        cov=load(name);validate_modes(cov)
        row={'rp_over_M':rp,'alpha':.3,'a_over_M':.8771530275949366}
        for b in ('infinity','horizon'):
            total=sum(v['flux'][b]['orbital_energy'] for v in cov[b]['modes'])
            l5=sum(v['flux'][b]['orbital_energy'] for v in cov[b]['modes'] if v['ell']<=5)
            ref=next(v for v in comparison['rows'] if v['rp']==rp and v['boundary']==b)
            np.testing.assert_allclose(abs(total),ref['local_finite_magnitude'],rtol=1e-10)
            np.testing.assert_allclose(abs(l5),ref['local_same_ellmax5_magnitude'],rtol=1e-10)
            row[b]={'signed_local_Dyson_cutoff_per_q2_eta':total,'signed_local_ell_le5_per_q2_eta':l5,'dyson_later_magnitude_per_q2_eta':ref['dyson_in_later_comparison_magnitude'],'li_magnitude_per_q2_eta':ref['li_magnitude'],'vs_Dyson_later_percent':100*(abs(total)/ref['dyson_in_later_comparison_magnitude']-1),'vs_Li_same_ellmax5_percent':100*(abs(l5)/ref['li_magnitude']-1),'computed_channels':len(cov[b]['modes'])}
        rows.append(row)
    partial=load('figure2_rp15_partial_20260916.json')
    for channel in partial['channels']:
        path=DATA/channel['file'];actual=json.loads(path.read_text());assert actual['flux']==channel['flux'];INPUTS[str(path.relative_to(ROOT))]=digest(path)
    save('figure2_rp15_partial_values.json',partial)
    fig,axs=plt.subplots(2,2,figsize=(11.8,7.5),layout='constrained',gridspec_kw={'height_ratios':[2.8,1]})
    for j,b in enumerate(('infinity','horizon')):
        ax=axs[0,j]
        for src,label,color in [('dyson','Dyson data in later comparison','#D55E00'),('li','Li et al. (a/M = 0.88)','#009E73')]:
            xy=curves[src+'_'+b]; ax.semilogy(xy[:,0],xy[:,1],color=color,label=label,lw=1.5)
        x=[r['rp_over_M'] for r in rows]
        ax.semilogy(x,[abs(r[b]['signed_local_Dyson_cutoff_per_q2_eta'])/.3**6 for r in rows],'o',color='#0072B2',label='Local: Dyson mode cutoff',zorder=4)
        if b=='infinity':
            ax.semilogy(x,[abs(r[b]['signed_local_ell_le5_per_q2_eta'])/.3**6 for r in rows],'s',mfc='white',mec='#CC79A7',label=r'Local: $\ell\leq5$ (Li cutoff)',zorder=5)
        if b=='horizon':
            assert not partial['boundaries'][b]['is_total']
            ax.semilogy([15],[abs(partial['boundaries'][b]['partial_signed_flux'])/.3**6],'^',mfc='white',mec='#6A3D9A',ms=8,label='New r=15M: 5/18 H modes only',zorder=6)
        ax.set(title='Infinity' if b=='infinity' else 'Horizon magnitude (local signed flux < 0)',xlim=(5,50),ylabel=r'$|F^s|/(q^2\epsilon^2)$')
        ax.grid(alpha=.2);ax.legend(fontsize=8,loc='best')
        ax=axs[1,j]
        for key,label,color,marker in [('vs_Dyson_later_percent','vs. later Dyson','#D55E00','o'),('vs_Li_same_ellmax5_percent','vs. Li, same cutoff','#009E73','s')]:
            ax.plot(x,[r[b][key] for r in rows],marker,color=color,linestyle='none',label=label)
        ax.axhline(0,c='black',lw=.7);ax.set(xlabel=r'$r_p/M$',ylabel='Difference [%]',xlim=(5,50));ax.grid(alpha=.2)
    axs[1,0].legend(fontsize=8)
    fig.suptitle(r'Fig. 2 comparison: $\alpha=0.3$, three finite totals plus a new partial horizon sum'+'\n'+r'$\epsilon^2=\alpha^6 M_c/M$; curves are published references, markers are local calculations',fontsize=12)
    plot_save(fig,'figure2_flux_comparison')
    save('figure2_flux_values.json',{'normalization':'per q^2 eta, eta=M_c/M; signed local orbital-effective scalar flux','rows':rows,'status':'basic_agreement_at_sampled_orbits_not_full_radial_scan','reference_provenance':'li_reference_comparison_20260916.json; original vector PDF, not raw author tables','limitations':['Local a/M=.8771530275949366 vs Li .88; cloud imaginary-frequency convention not aligned.','Finite numerical settings; no certified all-mode numerical error bar.','Local infinity ell<=6 vs Dyson; ell<=5 vs Li. Horizon ell<=5 for both.','PDF interpolation has no rigorous error bound.']})
    save('figure2_published_reference_curves.json',{'status':'digitized_reference_only_not_local_simulation','normalization':'per q^2 epsilon^2','curves':{k:[{'rp_over_M':float(x),'flux_magnitude':float(y)} for x,y in xy] for k,xy in curves.items()}})
    return rows

def modal_figure():
    d=load('figure6_progress_L18.json'); assert d['computed_count']==d['required_count']==36 and not d['missing_modes']
    rows=d['modes']
    for r in rows:
        p=DATA/r['file'];a=json.loads(p.read_text());INPUTS[str(p.relative_to(ROOT))]=digest(p)
        np.testing.assert_allclose(a['flux']['infinity']['orbital_energy'],r['computed'],rtol=1e-14)
    fig,axs=plt.subplots(2,1,figsize=(10.5,7.1),layout='constrained',gridspec_kw={'height_ratios':[3,1]})
    colors=plt.get_cmap('turbo')
    for m in sorted({r['m'] for r in rows}):
        rs=[r for r in rows if r['m']==m];c=colors((m-2)/10)
        axs[0].semilogy([r['ell'] for r in rs],[r['computed']/.3**6 for r in rs],'o-',c=c,lw=.8,ms=4,label=f'm={m}')
        refs=[r for r in rs if r['digitized'] is not None]
        axs[0].semilogy([r['ell'] for r in refs],[r['digitized']/.3**6 for r in refs],'D',mfc='none',mec=c,ms=6)
        axs[1].plot([r['ell'] for r in refs],[100*(r['ratio']-1) for r in refs],'o',c=c,ms=4)
    axs[0].set(ylabel=r'$F^{s,\infty}_{\ell m}/(q^2\epsilon^2)$',ylim=(1e-21,.06),xlim=(1.7,12.3),xticks=range(2,13))
    axs[0].legend(ncol=6,fontsize=8,loc='lower left');axs[0].grid(alpha=.2)
    axs[1].set(xlabel=r'$\ell$',ylabel='Local / old ref. - 1 [%]',xlim=(1.7,12.3),xticks=range(2,13));axs[1].axhline(0,c='black',lw=.7);axs[1].grid(alpha=.2)
    fig.suptitle(r'Fig. 6: $r_p=20M$, 36 computed radiative modes through $\ell=12$'+'\nFilled circles: local; open diamonds: original Dyson markers (no corrected modal table available)',fontsize=12)
    plot_save(fig,'figure6_modal_flux')
    cumulative=[]
    for l in range(2,13):
        v=sum(r['computed'] for r in rows if r['ell']<=l)
        cumulative.append({'ellmax':l,'finite_infinity_flux_per_q2_eta':v})
    save('figure6_modal_values.json',{'status':'complete_finite_mode_range_not_full_convergence','normalization':'per q^2 eta','modes':rows,'cumulative':cumulative,'limitations':['Original mode reference has not been replaced by corrected author modal data.','Agreement of the total with Li does not establish mode-by-mode agreement.','No arbitrary 5% renormalization applied.']})
    return cumulative

def field_figure():
    f1=load('wake_diagnostic_rp3.5_L18_sl12_mixed_wide_linear_eps.json')
    f5=load('figure5_threshold_wakes.json')
    assert not f1['missing_modes'] and len(f1['included_modes'])==88
    assert all(not d['missing_modes'] and len(d['included_modes'])==88 for d in f5['datasets'])
    # Verify the saved reconstruction inputs against their creation-time manifests.
    for meta in [f1]+f5['datasets']:
        for ref in meta['inputs']:
            p=ROOT/ref['file'] if '/' in ref['file'] else DATA/ref['file']
            assert digest(p)==ref['sha256'], f'Stale field input: {p}'
            INPUTS[str(p.relative_to(ROOT))]=digest(p)
    z1=array('wake_diagnostic_rp3.5_L18_sl12_mixed_wide_linear_eps.npz');z5=array('figure5_threshold_wakes.npz')
    re=np.geomspace(*f1['plot']['radial_edges'],f1['plot']['radial_cells']+1)
    ae=np.linspace(-np.pi,np.pi,f1['plot']['angular_cells']+1)
    np.testing.assert_allclose(np.sqrt(re[:-1]*re[1:]),z1['radial_centers'],rtol=1e-14)
    np.testing.assert_allclose((ae[:-1]+ae[1:])/2,z1['angle_centers'],atol=1e-15)
    a=f1['parameters']['a'];rh=1+np.sqrt(1-a*a)
    fig,axs=plt.subplots(2,2,figsize=(12,10.4),layout='constrained')
    topmax=max(float(np.max(abs(z1[k]))) for k in ['equatorial','meridional'])
    bmax=max(float(np.max(abs(z5[k]))) for k in ['field_416','field_418'])
    plots=[(axs[0,0],re[:,None]*np.cos(ae),re[:,None]*np.sin(ae),z1['equatorial'],topmax,135,'Fig. 1: equatorial, r_p=3.5M','Y/M',3.5),
           (axs[0,1],re[:,None]*np.sin(ae),re[:,None]*np.cos(ae),z1['meridional'],topmax,135,'Fig. 1: meridional, r_p=3.5M','Z/M',3.5)]
    for ax,s,rp in [(axs[1,0],'416',41.6),(axs[1,1],'418',41.8)]:
        r=z5['radial_edges_'+s][:,None];ang=z5['angular_edges'][None,:]
        plots.append((ax,r*np.cos(ang),r*np.sin(ang),z5['field_'+s],bmax,225,f'Fig. 5: equatorial, r_p={rp}M','Y/M',rp))
    stats=[]
    for ax,x,y,z,vmax,extent,title,yl,rp in plots:
        assert np.isfinite(z).all()
        im=ax.pcolormesh(x,y,abs(z),cmap='magma',norm=Normalize(0,vmax),shading='flat',rasterized=True)
        ax.add_patch(Circle((0,0),rh,color='black'));ax.plot(rp,0,'+',color='cyan',ms=6)
        ax.set(xlim=(-extent,extent),ylim=(-extent,extent),aspect='equal',xlabel='X/M',ylabel=yl,title=title)
        fig.colorbar(im,ax=ax,shrink=.77,label=r'$|\delta\Phi|/(q\epsilon)$')
        stats.append({'panel':title,'rp_over_M':rp,'sampled_peak':float(np.max(abs(z))),'plot_grid_rms_not_volume_norm':float(np.sqrt(np.mean(abs(z)**2)))})
    fig.suptitle(r'Computed field slices: $\alpha=0.3$, $2\leq\ell\leq12$ (88 allowed modes per orbit)'+'\nFinite reconstructions; absolute color scales, no fit to the paper; BL-label spatial embedding',fontsize=12)
    plot_save(fig,'figures1_5_field_slices')
    for src,dst in [('wake_diagnostic_rp3.5_L18_sl12_mixed_wide_linear_eps.npz','figure1_complex_field.npz'),('figure5_threshold_wakes.npz','figure5_complex_fields.npz')]:shutil.copyfile(DATA/src,OUT/dst)
    save('figures1_5_field_values.json',{'normalization':'delta Phi/(q epsilon); epsilon=alpha^3 sqrt(M_c/M)','time':0,'coordinates':'X=r sin(theta) cos(phi), Y=r sin(theta) sin(phi), Z=r cos(theta), with BL labels; original paper embedding unconfirmed','mode_range':[2,12],'allowed_modes_per_orbit':88,'rows':stats,'status':'computed_field_patterns_not_quantitative_pixelwise_paper_validation','limitations':['Mixed discretizations and finite source domains retained.','No original complex field arrays are available for a pixelwise error test.','Magnitude plots are not total cloud density or r Phi.']})
    return stats

def physical_flux_figure(fluxrows):
    sweep=load('gravitational_flux_sweep_L12.json')
    for r in sweep['rows']:
        p=DATA/r['input_file'];assert digest(p)==r['sha256'];INPUTS[str(p.relative_to(ROOT))]=digest(p)
        raw=json.loads(p.read_text())
        for b in ('infinity','horizon'):
            np.testing.assert_allclose(sum(v[b] for v in raw['modes']),r['totals'][b],rtol=1e-14,atol=0)
    q=1e-6;eta=.1;rows=[]
    for fr in fluxrows:
        g=next(v for v in sweep['rows'] if v['rp']==fr['rp_over_M'])
        r={'rp_over_M':fr['rp_over_M']}
        for b in ('infinity','horizon'):
            fs=q*q*eta*fr[b]['signed_local_Dyson_cutoff_per_q2_eta']
            fg=q*q*g['totals'][b]
            r[b]={'scalar_effective_orbital_signed':fs,'GW_signed':fg,'scalar_over_GW':fs/fg}
        rows.append(r)
    fig,axs=plt.subplots(1,2,figsize=(11.6,4.7),layout='constrained')
    for b,color,marker,label in [('infinity','#0072B2','o','Infinity'),('horizon','#D55E00','s','Horizon magnitude')]:
        rs=[r for r in sweep['rows'] if 8<=r['rp']<=32]
        axs[0].semilogy([r['rp'] for r in rs],[q*q*abs(r['totals'][b]) for r in rs],color=color,label='GW: '+label)
        axs[0].semilogy([r['rp_over_M'] for r in rows],[abs(r[b]['scalar_effective_orbital_signed']) for r in rows],marker,mfc='white',mec=color,ms=7,label='Scalar: '+label)
        axs[1].semilogy([r['rp_over_M'] for r in rows],[r[b]['scalar_over_GW'] for r in rows],marker,color=color,ms=7,label=label)
    axs[0].set(xlabel=r'$r_p/M$',ylabel='Physical flux magnitude',xlim=(8,32),title='Local scalar and full finite GW sums')
    axs[1].axhline(1,color='gray',lw=.8,ls='--');axs[1].set(xlabel=r'$r_p/M$',ylabel=r'$F_s/F_{GW}$',xlim=(8,32),title='Ratios at computed scalar orbits')
    for ax in axs:ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle(r'Fig. 4 companion calculation: $q=10^{-6}$, $M_c/M=0.1$, $\alpha=0.3$'+'\nExplicit mass normalization; this does not reproduce the inconsistent scaling of the original plot',fontsize=11)
    plot_save(fig,'figure4_physical_flux_companion')
    save('figure4_physical_flux_values.json',{'status':'local_physical_normalization_companion_not_original_figure4_replication','q':q,'eta':eta,'alpha':.3,'epsilon_squared':.3**6*eta,'scalar_factor':q*q*eta,'GW_factor':q*q,'scalar_cutoffs':{'infinity':6,'horizon':5},'GW_ellmax':12,'GW_m_signs':'both, no second doubling','rows':rows,'limitations':['Scalar reference discrepancies and numerical limits remain.','Paper caption epsilon^2=0.1 alpha^3 differs from eta=.1; our parameters are explicit.','Earlier cross-figure audit found old scalar scaling and single +22 GW fingerprints; not corrected here by empirical factors.','No inspiral or floating-orbit conclusion inferred from these stationary finite sums.']})
    return rows


def main():
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.titlesize':11,'savefig.facecolor':'white'})
    rows=flux_figure();cum=modal_figure();fields=field_figure();physical=physical_flux_figure(rows)
    INPUTS[str(Path(__file__).relative_to(ROOT))]=digest(Path(__file__))
    save('bundle_manifest.json',{'status':'audited_finite_results_not_claimed_complete_paper_reproduction','inputs_sha256':INPUTS,'outputs_sha256':{p.name:digest(p) for p in OUT.iterdir() if p.is_file() and p.name!='bundle_manifest.json'},'new_physical_solve_in_this_script':False,'reference_fitted_scale':False})
    print(json.dumps({'flux_rows':rows,'cumulative':cum,'fields':fields,'input_files_checked':len(INPUTS)},ensure_ascii=False))
if __name__=='__main__':main()
