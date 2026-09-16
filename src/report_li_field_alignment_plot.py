"""Compare fresh finite-mode fields with the later Li field figure, without fitting.

The two stated mode selections in the reference are shown separately. Original
pixel values remain raster-derived estimates, not author numerical arrays.
"""
import argparse, csv, hashlib, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,ListedColormap
from PIL import Image
from environment_wake import EnvironmentalWake
from source_provenance import source_fingerprint,validate_saved_samples,local_dependency_hashes
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'docs/field_alignment_20260917'
D=ROOT/'docs/environment_reproduction'
POS={(2,2),(3,3),(4,2),(4,4),(5,3),(5,5)}
ALL={(l,m) for l in range(2,6) for m in range(-l,l+1) if (l+m)%2==0}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 pa=argparse.ArgumentParser();pa.add_argument('--extra',nargs='*',type=Path,default=[]);args=pa.parse_args()
 current=source_fingerprint();implementation=local_dependency_hashes(ROOT/'src',[Path(__file__).stem]);inputs={};selected={41.1:{},42.1:{}}
 OUT.mkdir(parents=True,exist_ok=True)
 paths=sorted((OUT/'nonstatic').glob('*.json'))+args.extra
 for path in paths:
  obj=json.loads(path.read_text());p=obj.get('parameters',{})
  if 'samples' not in obj or 'scalar_ell' not in p:continue
  orbit=p['metric']['orbital_radius'];key=(p['scalar_ell'],p['scalar_m'])
  if orbit not in selected or key not in ALL:continue
  if not np.isclose(p['alpha'],.3,rtol=0,atol=1e-14) or not np.isclose(p['metric']['a'],.8771530275949366,rtol=0,atol=1e-12):
   raise ValueError(f'Unexpected cloud parameters for this fixed comparison: {path}')
  if obj.get('status')!='truncated_single_mode_not_converged' or not obj.get('samples') or 'flux' not in obj:
   raise ValueError(f'Incomplete source response: {path}')
  validate_saved_samples(obj,current)
  if key in selected[orbit]:raise ValueError(f'Duplicate mode {orbit} {key}')
  selected[orbit][key]=obj;inputs[str(path.relative_to(ROOT))]=sha(path)
 for orbit,modes in selected.items():
  if set(modes)!=ALL:raise ValueError(f'Missing fresh modes at {orbit}: {sorted(ALL-set(modes))}')
 refpath=D/'li_field_pixel_reference_20260917.json';ref=json.loads(refpath.read_text());inputs[str(refpath.relative_to(ROOT))]=sha(refpath)
 original=ROOT/'outputs/paper_original_reference/li_2507_02045v2/source/scalar_radiation_11.png';rgb=np.asarray(Image.open(original));inputs[str(original.relative_to(ROOT))]=sha(original)
 if sha(original)!=ref['inputs_sha256'][str(original.relative_to(ROOT))]:raise ValueError('Reference PNG differs from the calibrated image')
 bar=ref['calibration'];bar_rgb=rgb[bar['bar_row_top']:bar['bar_row_bottom']+1,bar['bar_column'],:3][::-1]
 # The published colorbar is discrete and nonuniform in standard viridis.
 # Use its actual lookup table for a numerical color comparison.
 reference_cmap=ListedColormap(bar_rgb/255.,name='Li_original_colorbar')
 fig,axes=plt.subplots(2,3,figsize=(15,10.2),layout='constrained')
 points=[];stats=[];arrays={};ringrows=[];norm=Normalize(0,bar['nominal_colorbar_maximum'])
 polarpath=D/'li_field_polar_reference_20260917.npz';polar=np.load(polarpath);inputs[str(polarpath.relative_to(ROOT))]=sha(polarpath)
 if sha(polarpath)!=ref['polar_array']['sha256']:raise ValueError('Polar reference differs from the calibrated NPZ')
 expected_shape=(len(polar['rp']),len(polar['radii']),len(polar['phi']))
 for key in ['nominal_abs_field','lower_color_interval','upper_color_interval','RGB_distance']:
  if polar[key].shape!=expected_shape or not np.all(np.isfinite(polar[key])):raise ValueError(f'Invalid polar reference: {key}')
 if (not np.all(polar['lower_color_interval']<=polar['nominal_abs_field'])
     or not np.all(polar['nominal_abs_field']<=polar['upper_color_interval'])
     or np.any(polar['RGB_distance']>3)):raise ValueError('Unaccepted polar reference pixels')
 fig_ring,ax_ring=plt.subplots(2,3,figsize=(13,7.5),layout='constrained')
 for row,orbit in enumerate(selected):
  allwake=EnvironmentalWake(list(selected[orbit].values()),ellmax=5,allow_mixed_discretization=True)
  poswake=EnvironmentalWake([v for k,v in selected[orbit].items() if k in POS],ellmax=5,allow_partial=True,allow_mixed_discretization=True)
  a=allwake.parameters['a'];horizon=1+np.sqrt(1-a*a)
  edges=np.geomspace(horizon+5e-4,200,361);angles=np.linspace(-np.pi,np.pi,721)
  r=np.sqrt(edges[:-1]*edges[1:])[:,None];phi=(angles[:-1]+angles[1:])[None,:]/2
  panel=next(x for x in ref['calibration']['panels'] if x['rp']==orbit)
  cx,cy,rad=panel['cx'],panel['cy'],panel['radius_pixels']
  crop=rgb[int(round(cy-rad)):int(round(cy+rad)),int(round(cx-rad)):int(round(cx+rad))]
  axes[row,0].imshow(crop,extent=(-200,200,-200,200),origin='upper')
  axes[row,0].set_title(f'Li reference: orbit {orbit}M')
  orbitpoints=[p for p in ref['points'] if p['rp']==orbit]
  for col,(label,wake) in enumerate([('All 18 allowed modes',allwake),('Six positive m modes',poswake)],1):
   amplitude_unit=wake.parameters['alpha']**3*np.sqrt(wake.parameters['cloud_mass'])
   if not np.isfinite(amplitude_unit) or amplitude_unit<=0:raise ValueError('Invalid cloud amplitude normalization')
   field=wake.evaluate(r,np.pi/2,phi)/amplitude_unit
   im=axes[row,col].pcolormesh(edges[:,None]*np.cos(angles)[None,:],edges[:,None]*np.sin(angles)[None,:],abs(field),cmap=reference_cmap,norm=norm,shading='flat',rasterized=True)
   axes[row,col].plot([0,orbit],[0,0],'.',color='black',ms=3)
   axes[row,col].set_title(f'{label}; orbit {orbit}M')
   tag=str(orbit).replace('.','p')+'_'+('all' if col==1 else 'positive')
   arrays[tag]=field;arrays[tag+'_r_edges']=edges;arrays[tag+'_phi_edges']=angles
   stats.append(dict(orbit=orbit,selection=label,mode_count=len(wake.modes),modes=sorted(wake.modes),maximum=float(abs(field).max()),parameters=wake.parameters))
   xy=np.array([[q['x_over_M'],q['y_over_M']] for q in orbitpoints]);rr=np.hypot(xy[:,0],xy[:,1]);ph=np.arctan2(xy[:,1],xy[:,0]);values=wake.evaluate(rr,np.pi/2,ph)/amplitude_unit
   ringr=polar['radii'][:,None];ringphi=polar['phi'][None,:]
   ringfield=abs(wake.evaluate(ringr,np.pi/2,ringphi)/amplitude_unit)
   ringindex=int(np.flatnonzero(polar['rp']==orbit)[0]);refvals=polar['nominal_abs_field'][ringindex]
   low=polar['lower_color_interval'][ringindex];high=polar['upper_color_interval'][ringindex]
   arrays[tag+'_ring_amplitude']=ringfield
   for k,radius in enumerate(polar['radii']):
    v=ringfield[k];refv=refvals[k];ref_rms=float(np.sqrt(np.mean(refv**2)))
    ringrows.append(dict(orbit=orbit,selection=label,radius=float(radius),inside_color_interval_fraction=float(np.mean((v>=low[k])&(v<=high[k]))),local_angular_rms=float(np.sqrt(np.mean(v**2))),reference_nominal_angular_rms=ref_rms,reference_rms_interval=[float(np.sqrt(np.mean(low[k]**2))),float(np.sqrt(np.mean(high[k]**2)))],normalized_amplitude_rms_difference=float(np.sqrt(np.mean((v-refv)**2))/ref_rms)))
    if radius in [50,100,150]:
     panelcol=[50,100,150].index(radius);ar=ax_ring[row,panelcol]
     if col==1:ar.fill_between(polar['phi'],low[k],high[k],color='0.7',alpha=.7,label='Li color interval')
     ar.plot(polar['phi'],v,label=label,lw=1.1)
     ar.set(title=f'Orbit {orbit}M, extraction r={radius:g}M',xlabel='Azimuth phi [rad]',ylabel='Field amplitude')
   for q,z in zip(orbitpoints,values):
    v=float(abs(z));paper=q.get('nominal_field_if_max_0p37') if q.get('pixel_accepted') else None
    interval=q.get('interval_from_colorband_and_tick_rounding') if paper is not None else None
    points.append(dict(**q,selection=label,local_real=float(z.real),local_imag=float(z.imag),local_amplitude=v,local_over_reference=None if not paper else v/paper,relative_difference=None if not paper else v/paper-1,reference_interval_lower=None if interval is None else interval[0],reference_interval_upper=None if interval is None else interval[1],inside_reference_color_interval=None if interval is None else bool(interval[0]<=v<=interval[1])))
  for ax in axes[row]:ax.set(xlim=(-200,200),ylim=(-200,200),aspect='equal',xlabel='X/M',ylabel='Y/M')
 if source_fingerprint()!=current or local_dependency_hashes(ROOT/'src',[Path(__file__).stem])!=implementation:
  raise ValueError('Source or assembly implementation changed during evaluation')
 fig.colorbar(im,ax=axes[:,1:],label=r'$|\delta\Phi|/(q\alpha^3\sqrt{M_c/M})$',extend='max',shrink=.9)
 fig.suptitle('Same orbit radii, scalar ell=2..5, time=0; no amplitude or phase fit\nLocal stationary spin 0.877153; Li spin 0.88. Metric ell<=6: finite resolution, not a convergence certificate.',fontsize=12)
 for suffix in ['png','svg']:
  path=OUT/f'li_field_alignment.{suffix}';fig.savefig(path,dpi=160)
  if suffix=='svg':path.write_text('\n'.join(x.rstrip() for x in path.read_text().splitlines())+'\n')
 plt.close(fig)
 ax_ring[0,0].legend(fontsize=8)
 fig_ring.suptitle('Fixed radius and azimuth; no phase or amplitude fit. Gray bands are raster color intervals.')
 fig_ring.savefig(OUT/'li_field_ring_comparison.png',dpi=160);plt.close(fig_ring)
 arrays['comparison_radii']=polar['radii'];arrays['comparison_phi']=polar['phi']
 np.savez_compressed(OUT/'li_field_alignment.npz',**arrays)
 result=dict(status='fresh_finite_resolution_comparison_not_exact_parameter_match',source_provenance=current,inputs_sha256=inputs,statistics=stats,points=points,ring_comparisons=ringrows,normalization=ref['normalization'],amplitude_fitted=False,phase_fitted=False,reference_mode_scope=ref['mode_scope'],limitations=['Reference source pixels have contour/colorbar/coordinate quantization; do not interpret percent differences as precise author-data errors.','Li uses a=.88 and a quasi-bound cloud; local background is an exactly stationary threshold cloud with a=.8771530275949366.','Metric ell<=6, theta order12 and radial order8/h32 require controlled refinement before convergence claims.','Both full-18 and positive-six selections are reported; no choice is fitted to the image.'],implementation_sha256=sha(Path(__file__)),assembly_dependencies_sha256=implementation,display_colorbar=dict(method='Exact reversed original colorbar RGB rows; no fitted rescaling',row_top=bar['bar_row_top'],row_bottom=bar['bar_row_bottom'],column=bar['bar_column'],nominal_maximum=bar['nominal_colorbar_maximum']))
 (OUT/'li_field_alignment.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 with (OUT/'li_field_alignment_points.csv').open('w',newline='') as stream:
  keys=['rp','x_over_M','y_over_M','selection','nominal_field_if_max_0p37','local_amplitude','local_over_reference','relative_difference','reference_interval_lower','reference_interval_upper','inside_reference_color_interval']
  w=csv.DictWriter(stream,fieldnames=keys,extrasaction='ignore',lineterminator='\n');w.writeheader();w.writerows(points)
 print(json.dumps(stats,indent=2));print(OUT/'li_field_alignment.png')
if __name__=='__main__':main()
