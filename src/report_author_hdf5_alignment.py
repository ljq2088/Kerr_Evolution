"""Reconcile documented author HDF5 trace basis; reuse saved local metric values."""
import hashlib,json
from pathlib import Path
import numpy as np
from pybhpt.swsh import Yslm
from environment_angular_diagnostic import DenseRealHarmonic
from paper_author_hdf5 import COMPONENTS,load_mode
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"docs/environment_reproduction"
def enc(z):
    z=np.asarray(z,complex);return np.stack([z.real,z.imag],axis=-1).tolist()
def dec(z):
    z=np.asarray(z);return z[...,0]+1j*z[...,1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def mixing(m,q):
    x,w=np.polynomial.legendre.leggauss(q);t=np.arccos(x)
    c=.6*m/(8**1.5+.6);B=np.zeros((5,5),complex)
    for ell in range(m,5):
        S=DenseRealHarmonic(0,ell,m,c)(t)
        for j in range(m,5): B[j,ell]=2*np.pi*np.sum(w*Yslm(0,j,m,t).conj()*S)
    return B

def errors(value,ref):
    result=[]
    for k,name in enumerate(COMPONENTS):
        scale=np.linalg.norm(ref[:,k]);err=np.linalg.norm(value[:,k]-ref[:,k])
        result.append(dict(component=name,reference_l2=float(scale),difference_l2=float(err),relative_l2=float(err/scale) if scale>1e-20 else None))
    return result

def main():
    writer=ROOT/"outputs/paper_metric_reference/ConorDyson_KerrLorenzMSF/GenerationCodes/Numerics-Asymptotics/h1-Radial-gen/MetricReconstructRadiative.wl"
    reader=ROOT/"outputs/paper_metric_reference/ConorDyson_KerrLorenzMSF/GenerationCodes/Numerics-Asymptotics/h1-mmode-gen/h1FunctionalityAsymp.wl"
    result=dict(status="completed_basis_aligned_saved_comparison",scope="Original author 2026 fork sample output; not 2025 environmental paper raw flux data; not an execution of author solver",basis=dict(first_nine="spin-weighted spherical",trace_raw="spin-0 spheroidal",trace_aligned="spin-0 spherical",formula="B[j,ell] = 2*pi * integral Y_0jm(theta,0)^* S_0ellm(theta;a*m*Omega) sin(theta) dtheta; q10_aligned=B@q10_raw",quadrature=[64,96],output_ell=list(range(5)),input_spheroidal_ell=list(range(5)),unfitted_scale_phase="No scale, phase, sign, or m substitution fitted"),evidence=dict(writer=str(writer.relative_to(ROOT)),writer_sha256=sha(writer),writer_trace_line=2144,writer_projectors_lines=[2133,2143],reader=str(reader.relative_to(ROOT)),reader_sha256=sha(reader),old_bin_trace_conversion_lines=[197,229],HDF5_reader_omits_trace_conversion_lines=[250,285]),rows=[],inputs=[])
    for m in [1,2,3]:
        p=OUT/f"author_hdf5_metric_m{m}_q24.json";old=json.loads(p.read_text());mode=load_mode(m)
        B=mixing(m,64);B96=mixing(m,96)
        result['inputs'].append(dict(m=m,path=str(p.relative_to(ROOT)),sha256=sha(p),raw_hdf5_path=str(mode.path.relative_to(ROOT)),raw_hdf5_sha256=sha(mode.path),mixing=enc(B),mixing_q64_q96_max_abs=float(np.max(abs(B-B96)))))
        for row in old['rows']:
            raw=dec(row['reference']);local=dec(row['local'])
            actual=mode.at_index(row['side'],row['reference_grid_index'])[1].T
            if not np.array_equal(actual,raw):raise ValueError('saved reference is not exact raw same-m HDF5')
            aligned=raw.copy();aligned[:,9]=B@raw[:,9]
            out=dict(m=m,r=row['r'],side=row['side'],reference_grid_index=row['reference_grid_index'],reference_raw=enc(raw),reference_aligned=enc(aligned),local=enc(local),errors_raw=errors(local,raw),errors_aligned=errors(local,aligned))
            if 'local_omitting_dipole_vector' in row:
                out['errors_aligned_omitting_local_dipole_vector_DIAGNOSTIC_ONLY']=errors(dec(row['local_omitting_dipole_vector']),aligned)
            result['rows'].append(out)
    p=OUT/'author_hdf5_metric_aligned_comparison.json';p.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for r in result['rows']:print(r['m'],r['side'],r['r'],[e['relative_l2'] for e in r['errors_aligned']])
    return result
if __name__=='__main__':main()
