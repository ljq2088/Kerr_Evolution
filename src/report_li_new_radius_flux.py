"""Match the later paper's ell<=5 infinity flux at the new field-plot radii."""
import json,hashlib
from pathlib import Path
import numpy as np
import pymupdf as fitz
from scipy.interpolate import PchipInterpolator
from report_li_reference_comparison import points
from source_provenance import source_fingerprint,validate_saved_samples
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'docs/field_alignment_20260917'
MODES=[(2,2),(3,3),(4,2),(4,4),(5,3),(5,5)]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 pdf=ROOT/'outputs/paper_original_reference/li_2507_02045v2/source/total_flux_11.pdf'
 drawings=fitz.open(pdf)[0].get_drawings()
 nx=np.polyfit([61.12561798095703,98.98464965820312,136.8436737060547,174.7027130126953,212.56173706054688],[10,20,30,40,50],1)
 ny=np.polyfit([148.2509307861328,110.11843872070312,71.98595428466797,33.85345458984375],[-4,-3,-2,-1],1)
 curves={}
 for idx,name in [(42,'Li'),(44,'Dyson_in_Li')]:
  xy=points(drawings[idx]);curve=np.column_stack([np.polyval(nx,xy[:,0]),np.polyval(ny,xy[:,1])]);assert np.all(np.diff(curve[:,0])>0);curves[name]=curve
 current=source_fingerprint();rows=[];inputs={str(pdf.relative_to(ROOT)):sha(pdf)}
 for orbit in [41.1,42.1]:
  modes=[]
  for l,m in MODES:
   tag=str(orbit).replace('.','p');path=OUT/'nonstatic'/f'rp{tag}_sl{l}_sm{m}_L6_q12_nr8_h32.json'
   d=json.loads(path.read_text());validate_saved_samples(d,current);inputs[str(path.relative_to(ROOT))]=sha(path)
   modes.append(dict(ell=l,m=m,omega=d['parameters']['omega'],orbital_energy_infinity=d['flux']['infinity']['orbital_energy']))
  value=sum(x['orbital_energy_infinity'] for x in modes);references={}
  for label,curve in curves.items():
   reference=float(10**np.interp(orbit,curve[:,0],curve[:,1])*.3**6)
   pchip=float(10**PchipInterpolator(curve[:,0],curve[:,1])(orbit)*.3**6)
   references[label]=dict(reference=reference,local_over_reference=value/reference,relative_difference_percent=100*(value/reference-1),log_pchip_vs_linear_fraction=pchip/reference-1)
  rows.append(dict(orbit=orbit,local_infinity_flux=value,modes=modes,references=references))
 out=dict(status='fresh_L6_metric_finite_scalar_ell5_flux_comparison',rows=rows,normalization='per q^2 eta, eta=Mc/M; reference plotted q^2 zeta^2 converted by alpha^6',source_provenance=current,inputs_sha256=inputs,implementation_sha256=sha(Path(__file__)),limitations=['Local stationary a=.877153, reference a=.88; not an identical-background run.','Metric ell<=6, theta12, radial8/h32 is a finite resolution, not a converged production claim.','Reference is a vector-plot polyline readout, not author raw arrays; interpolation method difference is a sensitivity indicator, not an error bar.','No horizon total is formed because these field mode sets omit scalar ell0/1.'])
 (OUT/'li_new_radius_flux.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
 print(json.dumps(rows,indent=2))
if __name__=='__main__':main()
