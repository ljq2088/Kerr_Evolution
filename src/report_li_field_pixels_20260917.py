"""Digitize Li et al. Fig.7 field PNG on its printed colorbar, without local fitting."""
from pathlib import Path
import csv,hashlib,json
from PIL import Image
import numpy as np
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"docs/environment_reproduction"
PAPER=ROOT/"outputs/paper_original_reference/li_2507_02045v2/source"
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def main():
 p=PAPER/"scalar_radiation_11.png";a=np.asarray(Image.open(p));yr=np.arange(125,989)
 rgb=a[yr,1760,:3].astype(float);u=(988-yr)/(988-125);tree=cKDTree(rgb)
 def lookup(colors):
  ds,ind=tree.query(colors);lo=[];hi=[]
  for c,di,ix in zip(colors,ds,ind):
   exact=np.flatnonzero(np.all(rgb==c,axis=1))
   allowed=exact if len(exact) else tree.query_ball_point(c,3.0)
   lo.append(float(u[allowed].min()) if len(allowed) else np.nan)
   hi.append(float(u[allowed].max()) if len(allowed) else np.nan)
  return ds,np.array(lo),np.array(hi)
 panels=[dict(rp=41.1,cx=461.5,cy=556.5,radius_pixels=402.),dict(rp=42.1,cx=1315.,cy=556.5,radius_pixels=401.75)]
 rows=[];rings=[]
 for panel in panels:
  for x,y in [(-100,0),(100,0),(0,100),(0,-100),(-50,50),(50,-50),(-150,0),(150,0)]:
   ix=int(round(panel['cx']+panel['radius_pixels']*x/200));iy=int(round(panel['cy']-panel['radius_pixels']*y/200))
   color=a[iy-1:iy+2,ix-1:ix+2,:3].reshape(-1,3).astype(float);ds,lo,hi=lookup(color);ok=bool(max(ds)<=3)
   unit=float(np.median((lo+hi)/2)) if ok else None
   rows.append(dict(rp=panel['rp'],x_over_M=x,y_over_M=y,pixel_accepted=ok,max_RGB_distance=float(max(ds)),
     fraction_of_colorbar_maximum=unit,nominal_field_if_max_0p37=unit*.37 if ok else None,
     interval_from_colorband_and_tick_rounding=[float(lo.min()*.36875),float(hi.max()*.375)] if ok else None))
  for radius in [10,20,50,100,150,180]:
   ang=np.linspace(0,2*np.pi,1440,endpoint=False);ix=np.rint(panel['cx']+panel['radius_pixels']*radius*np.cos(ang)/200).astype(int);iy=np.rint(panel['cy']-panel['radius_pixels']*radius*np.sin(ang)/200).astype(int)
   ds,lo,hi=lookup(a[iy,ix,:3].astype(float));ok=ds<=3;v=(lo[ok]+hi[ok])/2*.37
   rings.append(dict(rp=panel['rp'],r_over_M=radius,accepted_fraction=float(ok.mean()),nominal_colorbar_maximum=.37,
    field_min=float(v.min()),field_median=float(np.median(v)),field_max=float(v.max()),angular_amplitude_rms=float(np.sqrt(np.mean(v*v)))))
 radii=np.array([10.,20.,30.,40.,50.,70.,100.,150.,180.]);phis=np.linspace(0,2*np.pi,72,endpoint=False)
 estimates=[];lower=[];upper=[];distance=[]
 for panel in panels:
  ix=np.rint(panel['cx']+panel['radius_pixels']*radii[:,None]*np.cos(phis)/200).astype(int);iy=np.rint(panel['cy']-panel['radius_pixels']*radii[:,None]*np.sin(phis)/200).astype(int)
  ds,lo,hi=lookup(a[iy,ix,:3].reshape(-1,3).astype(float))
  estimates.append(((lo+hi)*.37/2).reshape(len(radii),len(phis)));lower.append((lo*.36875).reshape(len(radii),len(phis)));upper.append((hi*.375).reshape(len(radii),len(phis)));distance.append(ds.reshape(len(radii),len(phis)))
 np.savez_compressed(OUT/'li_field_polar_reference_20260917.npz',rp=np.array([41.1,42.1]),radii=radii,phi=phis,nominal_abs_field=np.array(estimates),lower_color_interval=np.array(lower),upper_color_interval=np.array(upper),RGB_distance=np.array(distance))
 with (OUT/'li_field_pixel_reference_20260917.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 d=dict(status='Li_original_field_pixel_reference_different_parameters_not_local_agreement',source_url='https://arxiv.org/html/2507.02045v2',source_tex_lines=[148,164,534,535,626],
  normalization=dict(epsilon_Li='q=m_p/M',zeta_Li='alpha^3 sqrt(Mc/M)',figure_quantity='abs(Phi^(1,1))=abs(delta Phi)/(q alpha^3 sqrt(Mc/M))',no_additional_alpha_factor=True),
  parameters=dict(a=.88,alpha=.3,orbits=[41.1,42.1],time=0,r_max=200,scalar_ell_min=2,scalar_ell_max=5),
  mode_scope=dict(caption='all (ell,m) with 2<=ell<=5',body='modes that contribute to infinity flux',interpretation='Likely allowed positive m>=2: (2,2),(3,3),(4,2),(4,4),(5,3),(5,5); retain both textual statements as a scope ambiguity.'),
  calibration=dict(image_shape=list(a.shape),bar_column=1760,bar_row_top=125,bar_row_bottom=988,nominal_colorbar_maximum=.37,
   maximum_interval_from_rounded_tick_labels=[.36875,.375],printed_ticks=[.00,.07,.15,.22,.30,.37],panels=panels),
  points=rows,ring_statistics=rings,
  polar_array=dict(path='docs/environment_reproduction/li_field_polar_reference_20260917.npz',sha256=sha(OUT/'li_field_polar_reference_20260917.npz'),shape=[2,9,72],radii=radii.tolist(),phi_convention='0=right, pi/2=up, radians, no rotation or phase fit'),
  inputs_sha256={str(q.relative_to(ROOT)):sha(q) for q in [p,PAPER/'main.tex',Path(__file__)]},
  limitations=['Printed colorbar has rounded tick values; nominal values use .37 and carry <=1.4% endpoint uncertainty from labels.',
   'The image and colorbar are discretized into color bands. Exact matching RGB rows give a value interval, combined with rounded tick endpoint bounds; midpoints are nominal estimates only.',
   'Pixel inversion is not original author numerical data, nor proof that the plotted maximum equals the continuous field maximum.',
   'No local comparison performed because spin, orbits, and scalar cutoff differ.',
   'Circle centres and radii come from colored-pixel bounds; orientation is assumed as conventional x right and y up.',
   'Caption says equatorial theta=0 despite polar angular conventions. This is inconsistent wording, not permission to evaluate at the pole.',
   'Ring RMS refers to the magnitude displayed in the PNG, not the missing complex field phase.'])
 (OUT/'li_field_pixel_reference_20260917.json').write_text(json.dumps(d,indent=2,allow_nan=False)+'\n');print(json.dumps({'points':rows,'rings':rings},indent=2))

if __name__=='__main__':main()
