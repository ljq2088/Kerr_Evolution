"""Render coherent, complete field modes with the original fixed Fig.1 legend."""
import argparse,csv,json,hashlib
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicHermiteSpline
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,Normalize
from matplotlib.patches import Circle
from PIL import Image
from environment_source import angular_mode
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def required(lmax):return {(l,m) for l in range(2,lmax+1) for m in range(-l,l+1) if (l+m)%2==0}

def evaluate(modes,r,theta,phi,lmax):
    r,theta,phi=np.broadcast_arrays(r,theta,phi)
    total=np.zeros(r.shape,complex)
    for (ell,m),(omega,radial) in modes.items():
        if ell>lmax:continue
        with np.errstate(divide='ignore',invalid='ignore'):
            S=angular_mode(theta,ell,m,.88**2*(omega**2-.3**2))[0]
        total+=radial(r)*S*np.exp(1j*m*phi)
    return total/.3**3

def main():
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);args=p.parse_args();out=args.directory
    execution=json.loads((out/'execution.json').read_text())
    modes={};inputs={}
    for file in out.glob('mode_l*_m*.json'):
        d=json.loads(file.read_text());key=(d['ell'],d['m']);rp=out/d['radial_file']
        if d['config']!=execution['config'] or d['implementation_sha256']!=execution['implementation_sha256']:
            raise ValueError('Mixed field parameters/implementation')
        if key in modes or sha(rp)!=d['radial_sha256']:raise ValueError('Duplicate mode or changed radial data')
        with np.load(rp,allow_pickle=False) as z:
            radial=CubicHermiteSpline(z['r'],z['field'],z['derivative'],extrapolate=False)
        modes[key]=d['omega'],radial
        inputs[file.name]=sha(file);inputs[rp.name]=sha(rp)
    if set(modes)!=required(12):raise ValueError('All 88 field modes required')
    refpath=ROOT/'docs/environment_reproduction/figure1_pixel_comparison_20260917.json'
    ref=json.loads(refpath.read_text());cal=ref['calibration']
    source=ROOT/'outputs/paper_field_digitization_20260916/figure1_xref6.png'
    image=np.asarray(Image.open(source))
    lutfile=ROOT/'docs/environment_reproduction/figure1_original_legend_lut_20260917.csv'
    lut=np.genfromtxt(lutfile,delimiter=',',names=True)
    order=np.argsort(lut['field_value']);colors=np.column_stack([lut[c] for c in ('red','green','blue')])[order]/255
    cmap=ListedColormap(colors);norm=Normalize(0,cal['LUT_range'][1])
    # Check reference files against the earlier immutable calibration.
    for file in (source,lutfile):
        key=str(file.relative_to(ROOT))
        known={**ref.get('inputs_sha256',{}),**ref.get('outputs_sha256',{})}
        if key not in known or sha(file)!=known[key]:raise ValueError('Reference calibration changed')
    horizon=1+np.sqrt(1-.88**2)
    re=np.unique(np.r_[np.geomspace(horizon+.05,20,181),np.linspace(20,192,1501)])
    ae=np.linspace(-np.pi,np.pi,721);r=np.sqrt(re[:-1]*re[1:])[:,None];angle=((ae[:-1]+ae[1:])/2)[None,:]
    # Angular functions depend only on angle; evaluate them on the 1-D angular grid.
    fields={}
    for cutoff in (5,8,12):
        eq=np.zeros((len(r),angle.size),complex);mer=eq.copy()
        for (ell,m),(omega,radial) in modes.items():
            if ell>cutoff:continue
            with np.errstate(divide='ignore',invalid='ignore'):
                s_eq=angular_mode(np.pi/2,ell,m,.88**2*(omega**2-.3**2))[0]
                s_mer=angular_mode(np.arccos(np.cos(angle)),ell,m,.88**2*(omega**2-.3**2))[0]
            R=radial(r)/.3**3
            eq+=R*s_eq*np.exp(1j*m*angle)
            mer+=R*s_mer*np.exp(1j*m*np.where(np.sin(angle)>=0,0,np.pi))
        fields[cutoff]=(eq,mer)
    np.savez_compressed(out/'figure1_fields.npz',r=r[:,0],angle=angle[0],
        **{f'{plane}_l{cut}':arr for cut,pair in fields.items() for plane,arr in zip(('equatorial','meridional'),pair)})
    boxes={'equatorial':(519,2399,469,2350),'meridional':(519,2399,2391,3386)}
    fig,axes=plt.subplots(2,2,figsize=(12,10),layout='constrained',gridspec_kw={'height_ratios':[1,.53]})
    for i,key in enumerate(boxes):
        x0,x1,y0,y1=boxes[key]
        extent=[float(np.polyval(cal['x_polynomial'],x0-.5)),float(np.polyval(cal['x_polynomial'],x1-.5)),
                float(np.polyval(cal['y_polynomials'][key],y1-.5)),float(np.polyval(cal['y_polynomials'][key],y0-.5))]
        axes[i,0].imshow(image[y0:y1,x0:x1,:3],extent=extent,origin='upper')
        xx,yy=(re[:,None]*np.cos(ae),re[:,None]*np.sin(ae)) if i==0 else (re[:,None]*np.sin(ae),re[:,None]*np.cos(ae))
        im=axes[i,1].pcolormesh(xx,yy,abs(fields[12][i]),cmap=cmap,norm=norm,shading='flat',rasterized=True)
        axes[i,1].add_patch(Circle((0,0),horizon,color='black'))
        axes[i,1].plot(3.5,0,'k.',ms=3)
        axes[i,0].set_title('Dyson Fig.1: '+key)
        axes[i,1].set_title('Fresh a=0.88, scalar ell=2..12: '+key)
        for ax in axes[i]:ax.set(xlim=extent[:2],ylim=extent[2:],aspect='equal',xlabel='X/M',ylabel='Y/M' if i==0 else 'Z/M')
    fig.colorbar(im,ax=axes[:,1],label='Absolute scalar perturbation / (q alpha^3 sqrt(Mc/M))')
    fig.savefig(out/'figure1_comparison.png',dpi=170);plt.close(fig)
    rows=[]
    for row in ref['points']:
        x,y=row['x_over_M'],row['y_or_z_over_M'];r0=np.hypot(x,y)
        if row['plane']=='equatorial':theta=np.pi/2;phi=np.arctan2(y,x)
        else:theta=np.arccos(y/r0);phi=0 if x>=0 else np.pi
        values={str(l):float(abs(evaluate(modes,r0,theta,phi,l))) for l in (5,8,12)}
        rows.append(dict(plane=row['plane'],x=x,y_or_z=y,local=values,
            paper=row['paper_value'],paper_interval=row['paper_RGB_tolerance_interval'],
            ratio=values['12']/row['paper_value'] if row['pixel_accepted'] and row['paper_value'] else None))
    convergence={}
    for cut in (5,8):
        convergence[str(cut)]={key:float(np.linalg.norm(fields[cut][i]-fields[12][i])/np.linalg.norm(fields[12][i]))
            for i,key in enumerate(('equatorial','meridional'))}
    save(out/'figure1_comparison.json',dict(status='finite_resolution_field_comparison_not_convergence_claim',
        points=rows,scalar_truncation_complex_grid_l2=convergence,inputs_sha256=inputs,
        reference_calibration_sha256=sha(refpath),reference_image_sha256=sha(source),
        conventions=execution['config'],no_amplitude_or_phase_fit=True,
        limitations=execution['limitations']+['Pixel colors are not author raw arrays; paper coordinate embedding remains unverified.',
            'Truncation L2 compares the stored polar grid without physical volume weights.']))
    # Copy final human-facing outputs into the task workspace.
    target=Path('/mnt/c/Users/赖景祺/Documents/ChatGPT/EMRI环境效应/FIG1复现结果')
    target.mkdir(parents=True,exist_ok=True)
    import shutil
    for name in ('figure1_comparison.png','figure1_comparison.json','figure1_fields.npz'):
        shutil.copy2(out/name,target/name)
if __name__=='__main__':main()
