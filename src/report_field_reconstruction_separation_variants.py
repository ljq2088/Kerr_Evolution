"""Controlled wrong-separation-constant variants on a frozen Fig5 scalar22 source.

No variant is a production fix or an inferred historical author bug. Only the
radial A value and consistently its homogeneous boundary expansions change.
Source, angular projection, field angular factor, frequency and flux formula
remain frozen. This quantifies sensitivity to the published error category.
"""
import hashlib,json,time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import environment_radial as radial
from environment_response import SampledResponse
from environment_source import angular_mode
from environment_cloud import mode_flux
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'docs/environment_reproduction'
def enc(z):
    a=np.asarray(z,complex);return np.stack([a.real,a.imag],axis=-1).tolist()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    output=D/'field_reconstruction_separation_variants_20260917.json';ref=json.loads((D/'figure5_threshold_wakes.json').read_text());original=radial.angular_eigenvalue;result=dict(status='in_progress_frozen_source_wrong_equation_control',inputs={},rows=[],equation='(Delta Rprime)prime + [K^2/Delta - mu^2*r^2 - A - (a^2*omega^2-2*a*m*omega) - dA]R = J',limitations=['Frozen source and angular function make each nonzero dA an intentionally inconsistent equation experiment, not another mathematically equivalent convention.','The later paper states an original separation-constant error but does not specify its exact old expression. No tested variant is selected as the historical truth.','No production files, source data, waveform phase or flux prefactors are changed.','Historical Fig5 source reports do not carry source-code provenance; input bytes are hashed.']);allcurves=[];start=time.perf_counter()
    def save():output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for dataset in ref['datasets']:
        for item in dataset['inputs']:
            path=ROOT/item['file'];data=json.loads(path.read_text());p=data['parameters']
            if (p['scalar_ell'],p['scalar_m'])==(2,2):break
        result['inputs'][str(path.relative_to(ROOT))]=sha(path);a=p['metric']['a'];w=p['omega'];mu=p['alpha'];m=2;omega_c=w-1/(p['metric']['orbital_radius']**1.5+a);shift=a*a*w*w-2*a*m*w;A=original(2,2,a*a*(w*w-mu*mu));S=angular_mode(np.pi/2,2,2,a*a*(w*w-mu*mu))[0]
        variants=[('correct',0.,'correct A'),('double_lambda_shift',shift,'A -> A + (a^2*w^2-2amw): double-count affine shift'),('omit_lambda_shift',-shift,'A -> A - (a^2*w^2-2amw): remove affine shift'),('plus_a2mu2',a*a*mu*mu,'A -> A+a^2*mu^2'),('minus_a2mu2',-a*a*mu*mu,'A -> A-a^2*mu^2'),('plus_a2omega2',a*a*w*w,'A -> A+a^2*w^2'),('minus_a2omega2',-a*a*w*w,'A -> A-a^2*w^2')]
        rows=[];rgrid=np.linspace(3,225,445);rpoint=np.array([20.,50.,100.,200.]);baseline=None;curves={}
        for name,delta,meaning in variants:
            radial.angular_eigenvalue=lambda *args,_delta=delta,**kw:original(*args,**kw)+_delta
            try:response=SampledResponse.from_report(data)
            finally:radial.angular_eigenvalue=original
            g=response.green;values=response.evaluate(rpoint)[0];curve=response.evaluate(rgrid)[0]*S/mu**3;curves[name]=curve;zin=response.up_coefficient if g.propagating else 0j;zh=response.horizon_coefficient;flux=mode_flux(w,m,omega_c,1,mu,a,zin,zh);peaks,properties=find_peaks(-abs(curve),prominence=.005*np.max(abs(curve)))
            # Saved Gauss weights provide a separate response-coefficient consistency indicator.
            rn=np.array([s['r'] for s in data['samples']]);wg=np.array([s['weight'] for s in data['samples']]);J=np.array([complex(*s['source']) for s in data['samples']]);saved_grid_h=np.dot(wg,g.upsol.sol(rn)[0]*J)/g.w0
            row=dict(name=name,meaning=meaning,dA=float(delta),effective_A=float(np.real(A+delta)),effective_radial_lambda=float(np.real(A+shift+delta)),radii=rpoint.tolist(),radial_field=enc(values),equatorial_component_per_qepsilon=enc(values*S/mu**3),ZH=enc(zh),Zinf=enc(zin),flux=flux,relative_continuous_vs_saved_gauss_ZH=abs(zh-saved_grid_h)/max(abs(saved_grid_h),1e-300),wronskian_relative_spread=float(np.max(abs(g.wronskian(rpoint)/g.w0-1))),radial_modulus_minima_r=rgrid[peaks].tolist(),field_curve=enc(curve))
            if baseline is None:baseline=row
            else:
                b=np.array([complex(*v) for v in baseline['radial_field']]);row['field_amplitude_ratio_to_correct']=(abs(values)/np.maximum(abs(b),1e-300)).tolist();row['complex_field_relative_difference']=(abs(values-b)/np.maximum(abs(b),1e-300)).tolist();row['horizon_flux_ratio_to_correct']=flux['horizon']['orbital_energy']/baseline['flux']['horizon']['orbital_energy'];row['infinity_flux_ratio_to_correct']=flux['infinity']['orbital_energy']/baseline['flux']['infinity']['orbital_energy'] if baseline['flux']['infinity']['orbital_energy'] else None
            rows.append(row);result['rows']=[r for r in result['rows'] if r['orbit']!=p['metric']['orbital_radius']]+[dict(orbit=p['metric']['orbital_radius'],omega=w,mu=mu,correct_angular_A=float(np.real(A)),affine_lambda_shift=shift,source_frozen=True,angular_frozen=True,radial_grid=rgrid.tolist(),variants=rows)];save();print(p['metric']['orbital_radius'],name,row.get('field_amplitude_ratio_to_correct'),row.get('horizon_flux_ratio_to_correct'),flush=True)
        allcurves.append((p['metric']['orbital_radius'],rgrid,curves))
    fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained');styles=[('correct','Correct A','black','-'),('double_lambda_shift','Double affine shift','#bd3b28','-'),('omit_lambda_shift','Omit affine shift','#287aba','-'),('plus_a2mu2','Add a^2 mu^2','#9467bd','--'),('minus_a2mu2','Subtract a^2 mu^2','#329b62','--')]
    for axrow,(r0,rgrid,curves) in zip(axes,allcurves):
        for name,label,color,ls in styles:
            axrow[0].plot(rgrid,abs(curves[name]),label=label,color=color,ls=ls,lw=1.2);axrow[1].plot(rgrid,curves[name].real,label=label,color=color,ls=ls,lw=1.2)
        for ax in axrow:ax.axvline(r0,color='.6',lw=.8);ax.grid(alpha=.2);ax.set_xlabel('r/M')
        axrow[0].set_title(f'rp={r0}: scalar22 modulus');axrow[1].set_title(f'rp={r0}: fixed-phase real component');axrow[0].set_ylabel('Equatorial mode / (q epsilon)');axrow[1].set_ylabel('Real equatorial mode / (q epsilon)')
    axes[0,0].legend(fontsize=8);fig.suptitle('Frozen source: deliberately wrong radial separation constants; no fit or production change');fig.savefig(D/'field_reconstruction_separation_variants_20260917.png',dpi=160);fig.savefig(D/'field_reconstruction_separation_variants_20260917.pdf');plt.close(fig)
    result.update(status='completed_frozen_source_wrong_equation_control',elapsed_seconds=time.perf_counter()-start,implementation_sha256=sha(Path(__file__)),production_hashes={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'src/environment_radial.py',ROOT/'src/environment_cloud.py',ROOT/'src/environment_response.py',ROOT/'src/environment_source.py']});save()
if __name__=='__main__':main()
