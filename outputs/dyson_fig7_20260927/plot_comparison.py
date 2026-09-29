"""Compare unrescaled particle multipoles with digitized Dyson Fig.7."""
import json,sys,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent

def plot(data_path,label,stem):
    data=json.loads(Path(data_path).read_text())
    marker_path=ROOT/'docs/environment_reproduction/paper_figure7_markers.json'
    paper=json.loads(marker_path.read_text())
    ref={p['ell']:p['plotted_value'] for p in paper['markers']}
    rows=[r for r in data['multipoles'] if r['complete']]
    if not rows:return
    ell=np.array([r['ell'] for r in rows])
    y=np.array([r['abs_particle_sum_per_epsilon_q'] for r in rows])
    radial=np.array([r['abs_radial_sum_per_epsilon_q'] for r in rows])
    x=np.array(sorted(ref));p=np.array([ref[i] for i in x])
    fig,axs=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    ax=axs[0]
    ax.loglog(x,p,'D',color='black',label='Dyson Fig.7 (digitized)')
    ax.loglog(ell,y,'o-',label='Computed coherent field at particle')
    ax.loglog(ell,radial,':',color='gray',label='Radial-only sum: definition control')
    fits={}
    for name,xx,yy,color in [('paper',x,p,'black'),('computed',ell,y,'C0')]:
        fits[name]={}
        for lower in [6,8]:
            mask=(xx>=lower)&(xx<=12)
            if mask.sum()<4:continue
            b,c=np.polyfit(np.log(xx[mask]),np.log(yy[mask]),1)
            fits[name][str(lower)]={'p':float(-b),'C':float(np.exp(c)),'ell_min':lower,'ell_max':12,'n':int(mask.sum())}
            if lower==8:
                grid=np.linspace(8,12,100)
                ax.loglog(grid,np.exp(c)*grid**b,'--',color=color,label=f'{name} fit: p={-b:.4f}')
    ax.set(xlabel=r'$\ell$',ylabel=r'$|B_\ell|/(\epsilon q)$',title=r'$r_p=20M,\ \alpha=0.3$; no fitted rescaling')
    ax.legend(fontsize=8);ax.grid(True,which='both',alpha=.2)
    ratio=y/np.array([ref[int(i)] for i in ell])
    axs[1].plot(ell,ratio,'o-');axs[1].axhline(1,color='black',ls='--')
    axs[1].set(xlabel=r'$\ell$',ylabel='Computed / paper',title='Absolute amplitude agreement');axs[1].grid(alpha=.2)
    fig.suptitle(label,fontsize=12)
    fig.savefig(OUT/(stem+'.png'),dpi=180);fig.savefig(OUT/(stem+'.pdf'));plt.close(fig)
    result={'label':label,'fits':fits,'comparison':[{'ell':int(i),'computed':float(a),'paper':ref[int(i)],'ratio':float(r),'relative_error':float(r-1)} for i,a,r in zip(ell,y,ratio)],'inputs_sha256':{str(p):hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in [data_path,marker_path]},'limitations':['Digitized markers are not author numerical data; no formal digitization uncertainty.','Fit coefficients are descriptive; no convergence or reproduction claim.','Physical angular-evaluated sum and literal radial sum are shown separately.']}
    (OUT/(stem+'.json')).write_text(json.dumps(result,indent=2))
    return result
if __name__=='__main__':
    plot(ROOT/'docs/environment_reproduction/particle_field_L18.json','HISTORICAL DATA DIAGNOSTIC — fresh calculation pending','historical_vs_paper')
