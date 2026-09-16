"""Fig.1 original-pixel comparison on the original fixed colorbar scale."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
from PIL import Image
import pymupdf as fitz
from scipy.spatial import cKDTree
from scipy.interpolate import RegularGridInterpolator
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.patches import Circle
from report_paper_field_pixels_20260917 import groups

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"docs/environment_reproduction"
RAW=ROOT/"outputs/paper_field_digitization_20260916"
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def main():
 pdf=ROOT/"outputs/paper_original_reference/arxiv_v1_source/Figures/BosonEMRIsFieldPlots.pdf"
 doc=fitz.open(pdf);RAW.mkdir(exist_ok=True)
 for entry in doc[0].get_images(full=True):
  base=fitz.Pixmap(doc,entry[0]);mask=fitz.Pixmap(doc,entry[1]);fitz.Pixmap(base,mask).save(RAW/f"figure1_xref{entry[0]}.png")
 a=np.asarray(Image.open(RAW/"figure1_xref6.png"));bar=np.asarray(Image.open(RAW/"figure1_xref7.png"))
 w=(bar[:,:,:3].min(axis=2)>245)&(bar[:,:,3]>250)
 ticks=[x for x in groups(np.flatnonzero(w[:,85:113].sum(axis=1)>20)) if len(x)>5]
 ticky=np.array([x.mean() for x in ticks]);tickv=np.array([.020,.015,.010,.005]);coef=np.polyfit(ticky,tickv,1)
 assert max(abs(np.polyval(coef,ticky)-tickv))<4e-6
 yl=np.arange(6,2514);values=np.polyval(coef,yl);rgb=bar[yl,35,:3].astype(float);tree=cKDTree(rgb);top=float(values.max())
 with (OUT/"figure1_original_legend_lut_20260917.csv").open("w",newline="") as f:
  ww=csv.writer(f);ww.writerow(["legend_row","field_value","red","green","blue"]);ww.writerows(zip(yl,values,*rgb.T))
 xt=np.array([750,1104.5,1458.5,1812.5,2167]);xcoef=np.polyfit(xt,[-100,-50,0,50,100],1)
 yt={"equatorial":[701,1055,1409,1763.5,2117.5],"meridional":[2533.5,2887.5,3242]}
 yv={"equatorial":[100,50,0,-50,-100],"meridional":[50,0,-50]}
 yc={k:np.polyfit(yt[k],yv[k],1) for k in yt}
 boxes={"equatorial":(519,2399,469,2350),"meridional":(519,2399,2391,3386)}
 npz=OUT/"wake_diagnostic_rp3.5_L18_sl12_mixed_wide_linear_eps.npz";data=np.load(npz);r=data["radial_centers"];angle=data["angle_centers"]
 ratio=r[1]/r[0];re=r[0]/np.sqrt(ratio)*ratio**np.arange(len(r)+1);ae=np.linspace(-np.pi,np.pi,len(angle)+1)
 cmap=LinearSegmentedColormap.from_list("fig1_legend",[(0,rgb[-1]/255)]+[(float(v/top),c/255) for v,c in zip(values[::-1],rgb[::-1])],N=2048)
 norm=Normalize(0,top);rows=[]
 fig,axes=plt.subplots(2,2,figsize=(11.2,9.6),layout="constrained",gridspec_kw={"height_ratios":[1,.53]})
 points={"equatorial":[(-100,0),(100,0),(0,100),(0,-100),(-50,50),(-50,-50),(50,-50),(20,20)],
         "meridional":[(-100,20),(100,20),(-50,20),(50,20),(-20,40),(20,40),(-50,-20),(50,-20)]}
 for rowid,key in enumerate(["equatorial","meridional"]):
  z=data[key];period=np.r_[angle[-1]-2*np.pi,angle,angle[0]+2*np.pi];zp=np.column_stack([z[:,-1],z,z[:,0]])
  inter=RegularGridInterpolator((np.log(r),period),zp,bounds_error=False,fill_value=np.nan)
  x0,x1,y0,y1=boxes[key];extent=[float(np.polyval(xcoef,x0-.5)),float(np.polyval(xcoef,x1-.5)),float(np.polyval(yc[key],y1-.5)),float(np.polyval(yc[key],y0-.5))]
  axes[rowid,0].imshow(a[y0:y1,x0:x1,:3],extent=extent,origin="upper",interpolation="nearest")
  axes[rowid,0].set_title("Original raster: "+key)
  xx,yy=(re[:,None]*np.cos(ae),re[:,None]*np.sin(ae)) if key=="equatorial" else (re[:,None]*np.sin(ae),re[:,None]*np.cos(ae))
  im=axes[rowid,1].pcolormesh(xx,yy,abs(z),norm=norm,cmap=cmap,shading="flat",rasterized=True)
  axes[rowid,1].add_patch(Circle((0,0),1.48021096,color="black"));axes[rowid,1].plot(3.5,0,"k.",ms=3)
  axes[rowid,1].set_title(f"Local, fixed legend; max={abs(z).max():.5f}")
  for ax in axes[rowid]:ax.set(xlim=extent[:2],ylim=extent[2:],aspect="equal",xlabel="X/M",ylabel="Y/M" if key=="equatorial" else "Z/M")
  for x,y in points[key]:
   px=(x-xcoef[1])/xcoef[0];py=(y-yc[key][1])/yc[key][0];ix,iy=int(round(px)),int(round(py))
   if key=="equatorial":assert not (ix>=1710 and iy<=1158),"Point in inset"
   patch=a[iy-2:iy+3,ix-2:ix+3,:3].astype(float);ds,ii=tree.query(patch.reshape(-1,3));valid=bool(max(ds)<=3)
   paper=float(np.median(values[ii])) if valid else None
   allowed=sorted({k for nn in tree.query_ball_point(patch.reshape(-1,3),3) for k in nn});iv=[float(values[allowed].min()),float(values[allowed].max())] if valid else None
   ang=np.arctan2(y,x) if key=="equatorial" else np.arctan2(x,y)
   local=float(abs(inter([[np.log(np.hypot(x,y)),ang]])[0]))
   rows.append(dict(plane=key,x_over_M=x,y_or_z_over_M=y,pixel_x=float(px),pixel_y=float(py),pixel_accepted=valid,
       max_RGB_distance=float(max(ds)),paper_value=paper,paper_RGB_tolerance_interval=iv,
       local_interpolated_abs_field=local,local_over_pixel=local/paper if paper else None))
 fig.colorbar(im,ax=axes[:,1],label=r"$|\phi^{(1,1)}|$: original fixed legend",shrink=.85)
 fig.suptitle("Fig. 1: original pixels and unchanged local fields",fontsize=14)
 fig.supxlabel("No amplitude or phase fit. Original inset is excluded from pixel comparisons; the paper's final right-hand X tick is mislabelled -100.",fontsize=9)
 image=OUT/"figure1_fixed_legend_20260917.png";fig.savefig(image,dpi=170);plt.close(fig)
 table=OUT/"figure1_pixel_comparison_20260917.csv"
 with table.open("w",newline="") as f:
  ww=csv.DictWriter(f,fieldnames=list(rows[0]));ww.writeheader();ww.writerows(rows)
 report=dict(status="fixed_original_legend_pixel_comparison_not_author_raw_data",calibration=dict(legend_tick_rows=ticky.tolist(),legend_tick_values=tickv.tolist(),legend_value_polynomial=coef.tolist(),LUT_range=[float(values.min()),top],x_tick_pixels=xt.tolist(),x_tick_values=[-100,-50,0,50,100],x_polynomial=xcoef.tolist(),y_tick_pixels=yt,y_tick_values=yv,y_polynomials={k:v.tolist() for k,v in yc.items()}),points=rows,
     inputs_sha256={str(p.relative_to(ROOT)):sha(p) for p in [npz,pdf,RAW/"figure1_xref6.png",RAW/"figure1_xref7.png",Path(__file__)]},
     outputs_sha256={str(p.relative_to(ROOT)):sha(p) for p in [image,table,OUT/"figure1_original_legend_lut_20260917.csv"]},
     limitations=["Printed-pixel inversion, not original field arrays; colors outside the observed legend rejected.","Original inset excluded. Original rightmost -100 label is visibly inconsistent with ordered ticks; physical +100 follows tick order and central origin, with no fit to local data.","Local Cartesian embedding as in cached reconstruction remains an assumption; original image does not define a full coordinate map.","Local complex polar values interpolated before modulus. No amplitude, phase, or angular mode fitting.","Colorbar endpoint is not a recovered global field maximum."])
 (OUT/"figure1_pixel_comparison_20260917.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
 print(json.dumps({"legend_range":report["calibration"]["LUT_range"],"points":rows},indent=2))

if __name__=="__main__":main()
