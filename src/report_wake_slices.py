"""Plot actual response slices with an explicit missing-mode manifest."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Circle
from environment_wake import EnvironmentalWake


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('coverage',type=Path)
    parser.add_argument('--ellmax',type=int,default=5)
    parser.add_argument('--allow-partial',action='store_true')
    parser.add_argument('--allow-mixed-discretization',action='store_true')
    parser.add_argument('--wide-linear',action='store_true',
        help='Wide stacked slices and central inset for comparison with Fig.1 geometry')
    parser.add_argument('--epsilon-q',action='store_true',
        help='Use epsilon=alpha^3 sqrt(Mc/M), with no amplitude fit')
    args=parser.parse_args()
    folder=args.coverage.parent;coverage=json.loads(args.coverage.read_text())
    sections=('field',) if 'field' in coverage else ('infinity','horizon')
    files=sorted({row['file'] for boundary in sections
                  for row in coverage[boundary]['modes'] if 'flux' in row and 2<=row['ell']<=args.ellmax})
    reports=[json.loads((folder/name).read_text()) for name in files]
    wake=EnvironmentalWake(reports,ellmax=args.ellmax,allow_partial=args.allow_partial,
                           allow_mixed_discretization=args.allow_mixed_discretization)
    rp=1+np.sqrt(1-wake.parameters['a']**2)
    inner=max(rp+.05,wake.radial_domain[0],
              max(row['numerical']['source_panels'][0] for row in wake.mode_provenance))
    requested_outer=np.hypot(135.,135.)*1.01 if args.wide_linear else 80.
    outer=min(requested_outer,wake.radial_domain[1])
    if args.wide_linear and outer<requested_outer:
        raise ValueError('Field domain does not cover the corners of the wide slices')
    if inner>=outer:raise ValueError('No radial range for the requested slices')
    re=np.geomspace(inner,outer,241 if args.wide_linear else 121)
    ae=np.linspace(-np.pi,np.pi,481 if args.wide_linear else 241)
    r=np.sqrt(re[:-1]*re[1:])[:,None];angle=(ae[:-1]+ae[1:])[None,:]/2
    equatorial=wake.evaluate(r,np.pi/2,angle)
    meridional=wake.evaluate(r,np.arccos(np.cos(angle)),np.where(np.sin(angle)>=0,0.,np.pi))
    scale=wake.parameters['alpha']**-3 if args.epsilon_q else 1.
    equatorial*=scale;meridional*=scale
    maximum=max(np.max(abs(equatorial)),np.max(abs(meridional)))
    if args.wide_linear:
        norm=Normalize(vmin=0,vmax=maximum)
        fig,axes=plt.subplots(2,1,figsize=(8,11),layout='constrained',
                              gridspec_kw={'height_ratios':[270,140]})
    else:
        norm=LogNorm(vmin=maximum*1e-5,vmax=maximum)
        fig,axes=plt.subplots(1,2,figsize=(11,5.7),layout='constrained')
    for ax,values,X,Y,label in (
        (axes[0],equatorial,re[:,None]*np.cos(ae),re[:,None]*np.sin(ae),'Equatorial: theta=pi/2'),
        (axes[1],meridional,re[:,None]*np.sin(ae),re[:,None]*np.cos(ae),'Meridional: y=0')):
        im=ax.pcolormesh(X,Y,abs(values),norm=norm,cmap='magma',shading='flat',rasterized=True)
        ax.add_patch(Circle((0,0),rp,color='black'))
        ax.plot(wake.parameters['r0'],0,'+',color='cyan',ms=8,mew=1.5)
        ax.set(aspect='equal',xlabel='x/M',ylabel='y/M' if ax is axes[0] else 'z/M',title=label)
        if args.wide_linear:
            ax.set_xlim(-135,135)
            ax.set_ylim((-135,135) if ax is axes[0] else (-70,70))
    if args.wide_linear:
        inset=axes[0].inset_axes([.63,.63,.35,.35])
        inset.pcolormesh(re[:,None]*np.cos(ae),re[:,None]*np.sin(ae),
            abs(equatorial),norm=norm,cmap='magma',shading='flat',rasterized=True)
        inset.add_patch(Circle((0,0),rp,color='black'))
        inset.plot(wake.parameters['r0'],0,'+',color='cyan',ms=6)
        inset.set(xlim=(-14,14),ylim=(-14,14),aspect='equal')
        inset.tick_params(labelsize=7,colors='white',direction='in',pad=-13)
        for spine in inset.spines.values():spine.set_color('white')
    color_label=r'$|\delta\Phi|/(q\epsilon)$' if args.epsilon_q else r'$|\delta\Phi|/[q\sqrt{M_c/M}]$'
    fig.colorbar(im,ax=axes,label=color_label)
    missing=', '.join(f'({ell},{m})' for ell,m in wake.missing) or 'none'
    title='Partial scalar response' if wake.missing else 'Finite-mode scalar response'
    if wake.mixed_discretization:title+=' (mixed discretization)'
    fig.suptitle(f"{title}: alpha={wake.parameters['alpha']}, r_p={wake.parameters['r0']}M, ell=2..{args.ellmax}\nMissing modes: {missing}",fontsize=12)
    cutoff=wake.parameters['metric_ellmax']
    suffix='_mixed' if wake.mixed_discretization else ''
    if args.wide_linear:suffix+='_wide_linear'
    if args.epsilon_q:suffix+='_eps'
    stem=folder/f"wake_diagnostic_rp{wake.parameters['r0']:g}_L{cutoff}_sl{args.ellmax}{suffix}"
    fig.savefig(stem.with_suffix('.png'),dpi=160,bbox_inches='tight');plt.close(fig)
    np.savez_compressed(stem.with_suffix('.npz'),radial_centers=r[:,0],angle_centers=angle[0],
                        equatorial=equatorial,meridional=meridional)
    report=dict(status=wake.status,parameters=wake.parameters,missing_modes=wake.missing,
        mixed_discretization=wake.mixed_discretization,mode_provenance=wake.mode_provenance,
        common_radial_domain=wake.radial_domain,
        included_modes=sorted(wake.modes),time=0,
        plot=dict(layout='wide_stacked_with_inset' if args.wide_linear else 'side_by_side',
            color_scale='linear' if args.wide_linear else 'logarithmic',
            normalization='per q epsilon' if args.epsilon_q else 'per q sqrt(Mc/M)',
            applied_scale=scale,amplitude_fitted=False,radial_edges=[inner,outer],
            radial_cells=len(re)-1,angular_cells=len(ae)-1,
            figure1_coordinate_mapping_confirmed=False),
        coordinates='BL-label embedding x=r sin(theta)cos(phi), y=r sin(theta)sin(phi), z=r cos(theta); not Kerr-Schild Cartesian coordinates',
        quantity=('Complex scalar response per q epsilon' if args.epsilon_q else
                  'Complex scalar response per q sqrt(Mc/M)')+'; not density, total cloud, or rPhi',
        inputs=[dict(file=name,sha256=hashlib.sha256((folder/name).read_bytes()).hexdigest()) for name in files],
        limitation='Each mode retains its recorded finite source grid. Missing modes are explicitly omitted. Mixed discretization requires separate convergence checks. Not a paper wake reproduction.')
    stem.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print(stem.with_suffix('.png'),wake.status,wake.missing,flush=True)


if __name__=='__main__':main()
