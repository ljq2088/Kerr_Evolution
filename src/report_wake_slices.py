"""Plot actual response slices with an explicit missing-mode manifest."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
from environment_wake import EnvironmentalWake


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('coverage',type=Path)
    parser.add_argument('--ellmax',type=int,default=5)
    parser.add_argument('--allow-partial',action='store_true')
    args=parser.parse_args()
    folder=args.coverage.parent;coverage=json.loads(args.coverage.read_text())
    files=sorted({row['file'] for boundary in ('infinity','horizon')
                  for row in coverage[boundary]['modes'] if 'flux' in row and 2<=row['ell']<=args.ellmax})
    reports=[json.loads((folder/name).read_text()) for name in files]
    wake=EnvironmentalWake(reports,ellmax=args.ellmax,allow_partial=args.allow_partial)
    rp=1+np.sqrt(1-wake.parameters['a']**2)
    re=np.geomspace(max(rp+.05,wake.parameters['source_panels'][0]),80.,121)
    ae=np.linspace(-np.pi,np.pi,241)
    r=np.sqrt(re[:-1]*re[1:])[:,None];angle=(ae[:-1]+ae[1:])[None,:]/2
    equatorial=wake.evaluate(r,np.pi/2,angle)
    meridional=wake.evaluate(r,np.arccos(np.cos(angle)),np.where(np.sin(angle)>=0,0.,np.pi))
    maximum=max(np.max(abs(equatorial)),np.max(abs(meridional)))
    norm=LogNorm(vmin=maximum*1e-5,vmax=maximum)
    fig,axes=plt.subplots(1,2,figsize=(11,5.7),layout='constrained')
    for ax,values,X,Y,label in (
        (axes[0],equatorial,re[:,None]*np.cos(ae),re[:,None]*np.sin(ae),'Equatorial: theta=pi/2'),
        (axes[1],meridional,re[:,None]*np.sin(ae),re[:,None]*np.cos(ae),'Meridional: y=0')):
        im=ax.pcolormesh(X,Y,abs(values),norm=norm,cmap='magma',shading='flat',rasterized=True)
        ax.add_patch(Circle((0,0),rp,color='black'))
        ax.plot(wake.parameters['r0'],0,'+',color='cyan',ms=8,mew=1.5)
        ax.set(aspect='equal',xlabel='x/M',ylabel='y/M' if ax is axes[0] else 'z/M',title=label)
    fig.colorbar(im,ax=axes,label=r'$|\delta\Phi|/[q\sqrt{M_c/M}]$')
    missing=', '.join(f'({ell},{m})' for ell,m in wake.missing) or 'none'
    title='Partial scalar response' if wake.missing else 'Finite-mode scalar response'
    fig.suptitle(f"{title}: alpha={wake.parameters['alpha']}, r_p={wake.parameters['r0']}M, ell=2..{args.ellmax}\nMissing modes: {missing}",fontsize=12)
    stem=folder/f"wake_diagnostic_rp{wake.parameters['r0']:g}_L{wake.parameters['metric_ellmax']}_sl{args.ellmax}"
    fig.savefig(stem.with_suffix('.png'),dpi=160);plt.close(fig)
    np.savez_compressed(stem.with_suffix('.npz'),radial_centers=r[:,0],angle_centers=angle[0],
                        equatorial=equatorial,meridional=meridional)
    report=dict(status=wake.status,parameters=wake.parameters,missing_modes=wake.missing,
        included_modes=sorted(wake.modes),time=0,
        coordinates='BL-label embedding x=r sin(theta)cos(phi), y=r sin(theta)sin(phi), z=r cos(theta); not Kerr-Schild Cartesian coordinates',
        quantity='Complex scalar response per q sqrt(Mc/M); not density, total cloud, or rPhi',
        inputs=[dict(file=name,sha256=hashlib.sha256((folder/name).read_bytes()).hexdigest()) for name in files],
        limitation='Same finite source grid as flux reports. Missing modes are explicitly omitted. Not a paper wake reproduction.')
    stem.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print(stem.with_suffix('.png'),wake.status,wake.missing,flush=True)


if __name__=='__main__':main()
