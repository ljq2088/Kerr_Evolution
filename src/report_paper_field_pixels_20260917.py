"""Read the original Fig.5 raster and legend; compare without amplitude fitting.

The scalar array in the PDF does not exist: the recovered values are explicitly
colorbar-inverted display pixels. Bright colors outside its LUT remain invalid.
"""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import pymupdf as fitz
from PIL import Image
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Circle

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"docs/environment_reproduction"
RAW=ROOT/"outputs/paper_field_digitization_20260916"
PDF=ROOT/"outputs/paper_original_reference/arxiv_v1_source/Figures/ResonanceBosonEMRIsFieldPlot.pdf"
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()


def groups(ids):
    return [x for x in np.split(ids,np.flatnonzero(np.diff(ids)>1)+1) if len(x)]


def main():
    RAW.mkdir(exist_ok=True)
    doc=fitz.open(PDF)
    for ref in (6,7):
        pix=fitz.Pixmap(doc,ref);mask=fitz.Pixmap(doc,doc[0].get_images(full=True)[ref-6][1])
        fitz.Pixmap(pix,mask).save(RAW/f"figure5_xref{ref}.png")
    field=np.asarray(Image.open(RAW/"figure5_xref6.png"))
    bar=np.asarray(Image.open(RAW/"figure5_xref7.png"))
    white=(bar[:,:,:3].min(axis=2)>245)&(bar[:,:,3]>250)
    tickgroups=[x for x in groups(np.flatnonzero(white[:,85:113].sum(axis=1)>20)) if len(x)>5]
    ytick=np.array([x.mean() for x in tickgroups]);ticks=.035-.005*np.arange(7)
    assert len(ytick)==len(ticks)==7
    coef=np.polyfit(ytick,ticks,1)
    assert np.max(abs(np.polyval(coef,ytick)-ticks))<1e-10
    ylut=np.arange(6,3521);vlut=np.polyval(coef,ylut);rgb=bar[ylut,35,:3].astype(float)
    tree=cKDTree(rgb)
    top=float(vlut.max());bottom=float(vlut.min())
    with (OUT/"figure5_original_legend_lut_20260917.csv").open("w",newline="") as f:
        w=csv.writer(f);w.writerow(["legend_row","field_value","red","green","blue"])
        w.writerows(zip(ylut,vlut,*rgb.T))
    # Tick-derived physical coordinates: no registration fitted to local data.
    xtick=np.array([705,1082,1458.5,1835,2212]);xvalue=np.array([-180,-90,0,90,180])
    xcoef=np.polyfit(xtick,xvalue,1)
    ycoords={"416":np.array([656,1032.5,1409,1786,2163]),
             "418":np.array([2575,2952,3329,3705,4082])}
    yvalues=np.array([180,90,0,-90,-180])
    coeffs={k:np.polyfit(v,yvalues,1) for k,v in ycoords.items()}
    box={"416":(519,2399,469,2350),"418":(519,2399,2389,4269)}
    saved=np.load(OUT/"figure5_threshold_wakes.npz")
    theta=(saved["angular_edges"][:-1]+saved["angular_edges"][1:])/2
    rows=[];stats=[]
    # Prespecified points avoid the secondary, ticks, and overlaid panel text.
    points=[(-100,0),(100,0),(0,100),(0,-100),(-50,50),(-50,-50),
            (50,50),(50,-50),(-150,0),(150,0),(0,150),(0,-150),
            (-20,20),(-20,-20),(20,20),(20,-20)]
    fig,axes=plt.subplots(2,2,figsize=(11,10.4),layout="constrained")
    vv=vlut[::-1];cc=rgb[::-1]/255
    nodes=[(0.,cc[0])]+[(float(v/top),c) for v,c in zip(vv,cc)]
    cmap=LinearSegmentedColormap.from_list("paper_observed_legend",nodes,N=2048)
    cmap.set_over(tuple(rgb[0]/255));norm=Normalize(0,top,clip=False)
    for rowid,key in enumerate(["416","418"]):
        data=saved[f"field_{key}"];edges=saved[f"radial_edges_{key}"]
        radii=np.sqrt(edges[:-1]*edges[1:])
        angle_periodic=np.r_[theta[-1]-2*np.pi,theta,theta[0]+2*np.pi]
        data_periodic=np.column_stack([data[:,-1],data,data[:,0]])
        interp=RegularGridInterpolator((np.log(radii),angle_periodic),data_periodic,bounds_error=False,fill_value=np.nan)
        x0,x1,y0,y1=box[key];crop=field[y0:y1,x0:x1,:3]
        # Audit each third opaque plot pixel. Exclude greys (text, tick, BH dots).
        sample=field[y0:y1:3,x0:x1:3]
        valid=(sample[:,:,3]>250)&(np.ptp(sample[:,:,:3].astype(float),axis=2)>30)
        pixel=sample[:,:,:3][valid].astype(float)
        dist,index=tree.query(pixel)
        good=dist<=3
        inferred=vlut[index]
        s=dict(rp=float(key)/10,full_cached_max=float(abs(data).max()),
               observed_legend_max=top,plot_sample_count=int(len(pixel)),
               fraction_with_RGB_distance_le_3=float(good.mean()),
               fraction_outside_observed_LUT=float((~good).mean()),
               recovered_pixel_quantiles={str(q):float(np.quantile(inferred[good],q)) for q in [.1,.5,.9,.99]},
               maximum_recoverable_in_legend_pixel=float(inferred[good].max()),
               maximum_is_not_recovered=True,
               local_polar_cells_above_legend_fraction=float((abs(data)>top).mean()))
        stats.append(s)
        ax=axes[rowid,0]
        extent=[float(np.polyval(xcoef,x0-.5)),float(np.polyval(xcoef,x1-.5)),
                float(np.polyval(coeffs[key],y1-.5)),float(np.polyval(coeffs[key],y0-.5))]
        ax.imshow(crop,extent=extent,origin="upper",interpolation="nearest")
        ax.set_title(f"Original PDF raster: $r_p={float(key)/10:.1f}M$")
        ax=axes[rowid,1]
        a=saved["angular_edges"][None,:];e=edges[:,None]
        im=ax.pcolormesh(e*np.cos(a),e*np.sin(a),abs(data),cmap=cmap,norm=norm,shading="flat",rasterized=True)
        rr=radii[:,None];aa=theta[None,:]
        ax.contour(rr*np.cos(aa),rr*np.sin(aa),abs(data),levels=[top],colors="#b91c1c",linewidths=.8)
        ax.add_patch(Circle((0,0),1.48021096,color="black"))
        ax.plot(float(key)/10,0,"k.",ms=3)
        ax.set_title(f"Local field, same fixed legend; max={abs(data).max():.4f}")
        for j,(x,y) in enumerate(points):
            px=float((x-xcoef[1])/xcoef[0]);py=float((y-coeffs[key][1])/coeffs[key][0])
            ix,iy=int(round(px)),int(round(py));patch=field[iy-2:iy+3,ix-2:ix+3,:3].astype(float)
            ds,ii=tree.query(patch.reshape(-1,3));values=vlut[ii]
            accepted=bool(np.max(ds)<=3)
            paper=float(np.median(values)) if accepted else None
            allowed=sorted({k for near in tree.query_ball_point(patch.reshape(-1,3),3.0) for k in near})
            interval=[float(vlut[allowed].min()),float(vlut[allowed].max())] if accepted else None
            z=interp([[np.log(np.hypot(x,y)),np.arctan2(y,x)]])[0]
            local=float(abs(z))
            rows.append(dict(rp=float(key)/10,x_over_M=x,y_over_M=y,
              pixel_x=px,pixel_y=py,paper_value=paper,
              paper_patch_min=float(values.min()) if accepted else None,
              paper_patch_max=float(values.max()) if accepted else None,
              max_RGB_distance=float(ds.max()),pixel_accepted=accepted,
              paper_RGB_tolerance_interval=interval,
              local_interpolated_abs_field=local,
              local_over_pixel=local/paper if accepted and paper>0 else None,
              amplitude_fitted=False))
        for ax in axes[rowid]:
            ax.set(xlim=(-225,225),ylim=(-225,225),aspect="equal",xlabel="X/M",ylabel="Y/M")
    fig.colorbar(im,ax=axes[:,1],label=r"$|\phi^{(1,1)}|$; fixed original legend calibration",extend="max",shrink=.82)
    fig.suptitle("Fig. 5: fixed color scale, unchanged local field values",fontsize=14)
    fig.supxlabel("Red contour: local field exceeds the original legend range. Original bright colors outside the LUT have no recovered maximum.",fontsize=9)
    figure=OUT/"figure5_fixed_legend_20260917.png";fig.savefig(figure,dpi=170);plt.close(fig)
    # Use the fixed correspondence at r=100: angular/radial interpolation never fits the image.
    with (OUT/"figure5_pixel_comparison_20260917.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    report=dict(status="fixed_original_legend_pixel_comparison_not_author_raw_data",
       inputs_sha256={str(p.relative_to(ROOT)):sha(p) for p in [PDF,OUT/"figure5_threshold_wakes.npz",Path(__file__)]},
       paper_pdf=dict(metadata=doc.metadata,embedded_file_count=doc.embfile_count(),xml_metadata=doc.get_xml_metadata(),
                      text_characters=len(doc[0].get_text()),vector_drawing_count=len(doc[0].get_drawings()),
                      raw_field_dimensions=list(field.shape),raw_legend_dimensions=list(bar.shape)),
       calibration=dict(legend_tick_rows=ytick.tolist(),legend_tick_values=ticks.tolist(),
          value_polynomial_slope_intercept=coef.tolist(),LUT_range=[bottom,top],LUT_RGB_distance_limit=3,
          x_tick_pixels=xtick.tolist(),x_tick_values=xvalue.tolist(),x_polynomial=xcoef.tolist(),
          y_tick_pixels={k:v.tolist() for k,v in ycoords.items()},y_tick_values=yvalues.tolist(),
          y_polynomials={k:v.tolist() for k,v in coeffs.items()},patch_side_pixels=5),
       panel_statistics=stats,points=rows,
       limitations=[
          "Recovering a value from the printed colorbar is digitization, not an original scalar array.",
          "Legend tick values were read from the visible labels; calibration uses only original pixels.",
          "Pixels more than 3 RGB units from the observed legend LUT are rejected, never extrapolated.",
          "The brightest field colors extend beyond the independent legend's LUT; its upper endpoint is not the global field maximum.",
          "Axes use original tick registration; BL-label embedding is assumed for local Cartesian comparison and remains unconfirmed in the paper.",
          "Local values are complex bilinear interpolations on the cached polar grid before modulus, no amplitude, phase, or mode-selection fitting.",
          "The fraction of local polar cells above the legend is not a Cartesian area fraction.",
          "Unequal definitions or source truncation can still affect pointwise values; this plot is not evidence of solver convergence."
       ],outputs_sha256={str(p.relative_to(ROOT)):sha(p) for p in [figure,OUT/"figure5_pixel_comparison_20260917.csv",OUT/"figure5_original_legend_lut_20260917.csv"]})
    (OUT/"figure5_pixel_comparison_20260917.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"calibration":report["calibration"],"statistics":stats,"points":rows},indent=2))

if __name__=="__main__":main()
