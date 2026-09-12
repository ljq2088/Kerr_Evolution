"""Compare complete finite-mode equatorial wakes across the mass threshold."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Circle
from environment_wake import EnvironmentalWake


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('coverage',nargs=2,type=Path)
    parser.add_argument('--ellmax',type=int,default=12)
    parser.add_argument('--radial-cells',type=int,default=320)
    parser.add_argument('--angular-cells',type=int,default=640)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.radial_cells<2 or args.angular_cells<4:
        raise ValueError('Insufficient plotting grid')
    datasets=[]
    for coverage_path in args.coverage:
        coverage=json.loads(coverage_path.read_text())
        files=sorted({row['file'] for row in coverage['field']['modes']
                      if 'flux' in row and 2<=row['ell']<=args.ellmax})
        paths=[coverage_path.parent/name for name in files]
        reports=[json.loads(path.read_text()) for path in paths]
        wake=EnvironmentalWake(reports,ellmax=args.ellmax,allow_mixed_discretization=True)
        if wake.parameters['cloud_mass']!=1:
            raise ValueError('This plot requires unit cloud mass normalization')
        horizon=1+np.sqrt(1-wake.parameters['a']**2)
        inner=max(horizon+.05,wake.radial_domain[0],
                  max(row['numerical']['source_panels'][0] for row in wake.mode_provenance))
        outer=np.hypot(225.,225.)*1.001
        source_outer=min(row['numerical']['source_panels'][-1] for row in wake.mode_provenance)
        if outer>min(wake.radial_domain[1],source_outer):
            raise ValueError('The complete square must fit inside the recorded source domain')
        edges=np.geomspace(inner,outer,args.radial_cells+1)
        angles=np.linspace(-np.pi,np.pi,args.angular_cells+1)
        radii=np.sqrt(edges[:-1]*edges[1:])[:,None]
        centers=(angles[:-1]+angles[1:])[None,:]/2
        values=wake.evaluate(radii,np.pi/2,centers)*wake.parameters['alpha']**-3
        if not np.all(np.isfinite(values)):
            raise ValueError('Nonfinite reconstructed field')
        datasets.append(dict(wake=wake,edges=edges,angles=angles,values=values,horizon=horizon,
            inputs=[dict(file=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()) for path in paths]))
    datasets.sort(key=lambda row:row['wake'].parameters['r0'])
    params=[row['wake'].parameters for row in datasets]
    if [p['r0'] for p in params]!=[41.6,41.8]:
        raise ValueError('Figure 5 requires r0=41.6 and 41.8')
    for key in ('alpha','a','cloud_mass','coordinates','flux_scaling'):
        if params[0][key]!=params[1][key]:raise ValueError('Different physics or normalization')
    maximum=max(float(np.max(abs(row['values']))) for row in datasets)
    if maximum<=0:raise ValueError('Vanishing field cannot define the color scale')
    fig,axes=plt.subplots(2,1,figsize=(7,12),layout='constrained')
    for ax,row in zip(axes,datasets):
        e=row['edges'][:,None];a=row['angles'][None,:]
        im=ax.pcolormesh(e*np.cos(a),e*np.sin(a),abs(row['values']),
                        norm=Normalize(0,maximum),cmap='viridis',shading='flat',rasterized=True)
        ax.add_patch(Circle((0,0),row['horizon'],color='black'))
        ax.plot(row['wake'].parameters['r0'],0,'+',color='black',ms=5)
        ax.set(xlim=(-225,225),ylim=(-225,225),aspect='equal',xlabel='X/M',ylabel='Y/M',
               title=f"Equatorial slice: r_p={row['wake'].parameters['r0']}M")
    fig.colorbar(im,ax=axes,label=r'$|\delta\Phi|/(q\epsilon)$')
    fig.suptitle('Threshold comparison: complete finite mode range, not converged')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(Path(str(args.output)+'.png'),dpi=160);plt.close(fig)
    np.savez_compressed(Path(str(args.output)+'.npz'),
        radial_edges_416=datasets[0]['edges'],radial_edges_418=datasets[1]['edges'],
        angular_edges=datasets[0]['angles'],field_416=datasets[0]['values'],field_418=datasets[1]['values'])
    report=dict(status='complete_finite_mode_range_not_paper_reproduction',ellmax=args.ellmax,
        coordinates='BL-label embedding X=r sin(theta) cos(phi), Y=r sin(theta) sin(phi); paper coordinate mapping unconfirmed',
        plane='theta=pi/2',time=0,normalization='per q epsilon, epsilon=alpha^3 sqrt(Mc/M)',
        amplitude_fitted=False,shared_linear_color_range=[0,maximum],extent=[-225,225],
        radial_cells=args.radial_cells,angular_cells=args.angular_cells,
        datasets=[dict(parameters=row['wake'].parameters,included_modes=sorted(row['wake'].modes),
                       missing_modes=row['wake'].missing,maximum=float(np.max(abs(row['values']))),
                       inputs=row['inputs']) for row in datasets],
        limitation='Finite grids, metric and scalar truncations, source cutoff and boundary construction need convergence checks. Shared color maximum is derived from computed data, not fitted to the paper.')
    Path(str(args.output)+'.json').write_text(json.dumps(report,indent=2)+'\n')
    print(Path(str(args.output)+'.png'),flush=True)


if __name__=='__main__':main()
