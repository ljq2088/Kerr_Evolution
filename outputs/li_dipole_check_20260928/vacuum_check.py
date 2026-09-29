import os,sys,json,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'));OUT=Path(__file__).resolve().parent
from environment_dense_metric import configure_metric_backend
configure_metric_backend(True)
from lorenz_metric import nonstatic_metric
from lorenz_tensor import lorenz_constraint,linearized_einstein
from source_provenance import local_dependency_hashes
os.nice(15)
rows=[]
for r,t in [(2.,1.1),(5.,1.1),(15.,1.1),(25.,1.1)]:
 start=time.time();g,h=nonstatic_metric(r,t,20.,.88,1,1,order=8)
 H=np.array([[z.value for z in row] for row in h]);C=np.array([z.value for z in lorenz_constraint(g,h)]);E=np.array([[z.value for z in row] for row in linearized_einstein(g,h)])
 # Dimensionless coordinate-scale diagnostics, not invariant error norms.
 scale=np.linalg.norm(H);row=dict(r=r,theta=t,h_norm=float(scale),lorenz_norm=float(np.linalg.norm(C)),einstein_norm=float(np.linalg.norm(E)),lorenz_scaled=float(r*np.linalg.norm(C)/scale),einstein_scaled=float(r*r*np.linalg.norm(E)/scale),elapsed=time.time()-start);rows.append(row);print(row,flush=True)
 (OUT/'vacuum_residual.json').write_text(json.dumps({'rows':rows,'parameters':{'a':.88,'r0':20.,'ell_g':1,'m_g':1,'jet_order':8},'limitations':['Vacuum field equations and gauge only, not particle jump normalization validation.','Coordinate-scaled norms are diagnostic, not invariant.'],'hashes':local_dependency_hashes(ROOT/'src',['lorenz_metric','environment_dense_metric'])},indent=2))
