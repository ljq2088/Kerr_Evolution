"""Schema, raw sample extraction and elementary QC of unmodified author HDF5."""
from pathlib import Path
import hashlib,json
import numpy as np
from paper_author_hdf5 import ROOT,load_mode,COMPONENTS,SPINS

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def enc(z):
    z=np.asarray(z,complex);return np.stack([z.real,z.imag],axis=-1).tolist()
def tortoise(r,a=.6):
    rp=1+np.sqrt(1-a*a);rm=1-np.sqrt(1-a*a)
    return r+(rp+rm)/(rp-rm)*(rp*np.log((r-rp)/2)-rm*np.log((r-rm)/2))
def main():
    ref=ROOT/'outputs/paper_metric_reference'
    manifest=json.loads((ref/'manifest.json').read_text())
    base=ref/'ConorDyson_KerrLorenzMSF/GenerationCodes/Numerics-Asymptotics'
    writer=base/'h1-Radial-gen/MetricReconstructRadiative.wl'
    reader=base/'h1-mmode-gen/h1FunctionalityAsymp.wl'
    run=base/'h1-Radial-gen/RunFile.wl'
    result=dict(evidence_class='original_author_dataset',is_original_author_run_by_us=False,is_raw_2025_environmental_flux_data=False,retrieval_manifest_sha256=sha(ref/'manifest.json'),repositories=[{k:r[k] for k in ['repository','commit']} for r in manifest['repositories']],generation_commit='Not embedded in files; downloaded repository commit is not proof of generation commit',files=[],schema=dict(root_attributes={},dataset_attributes={},complex_storage={'Re':'little-endian float64','Im':'little-endian float64'},shape_meaning=['output harmonic index ell=0..4','radial sample index'],component_order=list(COMPONENTS),spin_weights=list(SPINS),first_nine_basis='spin-weighted spherical',trace_basis='raw spin-0 spheroidal; see writer line2144 and old-bin reader line213',r_in_order='orbit outward endpoint first, then decreasing r',r_up_order='orbit inward endpoint first, then increasing r',branch_definition='In is sourced solution for r<=rp; Up is sourced solution for r>=rp, not raw homogeneous functions',fourier_mode='single complex positive m; no conjugate-m addition or factor2 applied'),parameters=dict(a=.6,rp=8.,central_mass_M=1,particle_mass_normalization='unit small-body mass in source; source jump at writer579 equals -16*pi*S/(u^t*Delta), no mass prefactor',omega_formula='m/(rp**1.5+a)',r_plus=1.8,radial_coordinate='Boyer-Lindquist r/M',grid='uniform tortoise rstar at step0.25M',lmax=4,lplot=4,n=4,angres=4,nterms=10,kapord=8,working_precision_config=100,accuracy_goal_config=20,stored_precision='machine64; ExportOutput uses N and packed arrays',xhor_config=.0001,rinf_config=4000.,rmax_config=1000.,horord_config=6,inford_config=7,config_caveat='configuration and filename corroborate parameters; no HDF5 metadata proves options at generation; fractional a in config is printed as multiline 3/-/5'),code_evidence=[dict(path=str(p.relative_to(ROOT)),sha256=sha(p)) for p in [writer,reader,run]],line_evidence=dict(BuildGrid=[419,492],source_and_mass=[526,580],spherical_projection=[2111,2144],ExportOutput=[2214,2229],old_bin_trace_transform=[197,229],HDF5_reader=[250,285],harmonic_synthesis=[348,393],tetrad_conversion=[401,456]),quality_warnings=['Do not copy LoadData.wl miLoad substitution m1 components8/9 from m2. Direct extraction always uses one actual m group.','Current fork generation copies use lmins2 for spin1 MST/grid in some locations, potentially dropping ell1,m1. HDF5 generation revision is unknown.','Stored trace is spheroidal; HDF5 reader does not perform conversion present in old bin reader.','These lmax4 sample files are not a converged high-ell reference and not the raw arrays used for 2025 flux figures.','Several first-nine projected components have sizeable orbit-side discontinuity. This can diagnose writer/truncation/path issues but is not by itself a unique cause.','Never interpret a common scalar magnitude across different tetrad components as a coordinate invariant norm.'])
    for m in [1,2,3]:
        mode=load_mode(m);entry=dict(m=m,path=str(mode.path.relative_to(ROOT)),sha256=sha(mode.path),bytes=mode.path.stat().st_size,ell=mode.ell.tolist(),omega=m/(8**1.5+.6),sides={},samples=[],orbit_continuity=[])
        for side in ['In','Up']:
            rr=mode.radii[side];cc=mode.coefficients[side];rs=tortoise(rr)
            entry['sides'][side]=dict(shape=list(cc.shape),r_first=float(rr[0]),r_last=float(rr[-1]),rstar_first=float(rs[0]),rstar_last=float(rs[-1]),rstar_step_min=float(np.min(np.diff(rs))),rstar_step_max=float(np.max(np.diff(rs))),finite=bool(np.isfinite(cc).all()),max_abs_by_component_ell=np.max(abs(cc),axis=2).tolist())
            for target in ([3.,6.,8.] if side=='In' else [8.,12.,20.]):
                index=int(np.argmin(abs(rr-target)));r,c=mode.at_index(side,index)
                entry['samples'].append(dict(side=side,target_radius=target,grid_index=index,radius=r,coefficients_component_ell_complex_pairs=enc(c)))
        a=mode.coefficients['In'][:,:,0];b=mode.coefficients['Up'][:,:,0]
        for k,name in enumerate(COMPONENTS):
            scale=max(np.max(abs(a[k])),np.max(abs(b[k])))
            entry['orbit_continuity'].append(dict(component=name,maximum_absolute_jump=float(np.max(abs(b[k]-a[k]))),scale_maximum_either_side=float(scale),relative_to_largest_coefficient=float(np.max(abs(b[k]-a[k]))/scale) if scale>1e-20 else None))
        result['files'].append(entry)
    out=ROOT/'docs/environment_reproduction/author_hdf5_schema.json';out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(out)
if __name__=='__main__':main()
