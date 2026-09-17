"""Single-channel counterfactual linked to fixed Li color intervals. No fitting."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from environment_wake import EnvironmentalWake
from environment_source import angular_mode
from environment_response import SampledResponse
from source_provenance import source_fingerprint,validate_saved_samples
BASE=Path(__file__).resolve().parent
OLD=ROOT/'docs/field_alignment_20260917'
D=ROOT/'docs/environment_reproduction'
POS={(2,2),(3,3),(4,2),(4,4),(5,3),(5,5)}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    fp=source_fingerprint();manifest=json.loads((OLD/'li_field_alignment.json').read_text())
    refpath=D/'li_field_polar_reference_20260917.npz';ref=np.load(refpath);calpath=D/'li_field_pixel_reference_20260917.json';cal=json.loads(calpath.read_text())
    assert sha(refpath)==cal['polar_array']['sha256']
    inputs={str(refpath.relative_to(ROOT)):sha(refpath),str(calpath.relative_to(ROOT)):sha(calpath)}
    selected={41.1:{},42.1:{}}
    for path,expected in manifest['inputs_sha256'].items():
        if not path.endswith('.json') or 'field_pixel_reference' in path:continue
        pth=ROOT/path;obj=json.loads(pth.read_text());p=obj.get('parameters',{})
        if 'scalar_ell' not in p or 'samples' not in obj:continue
        assert sha(pth)==expected
        validate_saved_samples(obj,fp)
        orbit=p['metric']['orbital_radius'];key=p['scalar_ell'],p['scalar_m']
        if key in selected[orbit]:raise ValueError('Duplicate scalar input')
        selected[orbit][key]=obj;inputs[path]=expected
    arrays={};rows=[];totals=[]
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,3,figsize=(13,7),layout='constrained')
    for i,orbit in enumerate([41.1,42.1]):
        assert len(selected[orbit])==18
        tag=str(orbit).replace('.','p');apath=BASE/f'li_appendix_omission_rp{tag}.json';ab=json.loads(apath.read_text());inputs[str(apath.relative_to(ROOT))]=sha(apath)
        r=ref['radii'][:,None];phi=ref['phi'][None,:];ival=int(np.flatnonzero(ref['rp']==orbit)[0])
        target=ref['nominal_abs_field'][ival];lo=ref['lower_color_interval'][ival];hi=ref['upper_color_interval'][ival]
        for selection,chosen in [('all18',selected[orbit]),('positive6',{k:v for k,v in selected[orbit].items() if k in POS})]:
            wake=EnvironmentalWake(list(chosen.values()),ellmax=5,allow_partial=selection=='positive6',allow_mixed_discretization=True)
            unit=wake.parameters['alpha']**3*np.sqrt(wake.parameters['cloud_mass'])
            field=wake.evaluate(r,np.pi/2,phi)/unit
            baseline22=wake.modes[(2,2)][1];g=baseline22.green;p=selected[orbit][(2,2)]['parameters']
            samples=ab['source_samples'];miss=SampledResponse(g,p['source_panels'],[x['r'] for x in samples],[complex(*x['omitted']) for x in samples],log_first=p['horizon_log_first_panel'])
            S=angular_mode(np.pi/2,2,2,wake.parameters['a']**2*(p['omega']**2-p['alpha']**2))[0]
            delta=miss.evaluate(ref['radii'])[0][:,None]*S*np.exp(2j*phi)/unit
            counter=field-delta
            arrays[f'{tag}_{selection}_correct']=field;arrays[f'{tag}_{selection}_counterfactual']=counter
            for j,rr in enumerate(ref['radii']):
                record=dict(orbit=orbit,selection=selection,radius=float(rr))
                for label,value in [('correct_direct',abs(field[j])),('counterfactual_only22',abs(counter[j]))]:
                    norm=np.sqrt(np.mean(target[j]**2))
                    record[label]=dict(D_RMS=float(np.sqrt(np.mean((value-target[j])**2))/norm),angular_rms=float(np.sqrt(np.mean(value**2))),rms_over_reference=float(np.sqrt(np.mean(value**2))/norm),inside_color_interval_fraction=float(np.mean((value>=lo[j])&(value<=hi[j]))))
                rows.append(record)
                if selection=='all18' and rr in [50,100,150]:
                    ax=axes[i,[50,100,150].index(rr)];ax.fill_between(ref['phi'],lo[j],hi[j],color='.8',label='Li color interval');ax.plot(ref['phi'],abs(field[j]),label='Direct source');ax.plot(ref['phi'],abs(counter[j]),ls='--',label='Only 22: printed omission');ax.set(title=f'Orbit {orbit} M, r={rr:g} M',xlabel='Azimuth [rad]',ylabel='Field amplitude')
            flux={side:sum(v['flux'][side]['orbital_energy'] for v in chosen.values()) for side in ['horizon','infinity']}
            newflux={side:flux[side]+ab['fluxes']['counterfactual_printed'][side]['orbital_energy']-ab['fluxes']['correct_direct'][side]['orbital_energy'] for side in flux}
            totals.append(dict(orbit=orbit,selection=selection,global_D_RMS_correct=float(np.linalg.norm(abs(field)-target)/np.linalg.norm(target)),global_D_RMS_counterfactual=float(np.linalg.norm(abs(counter)-target)/np.linalg.norm(target)),flux_correct=flux,flux_counterfactual_only22=newflux))
            print('COMPLETE',orbit,selection,totals[-1],flush=True)
    axes[0,0].legend(fontsize=8);fig.suptitle('Only scalar 22 source replaced by literal-Appendix counterfactual; all other modes fixed. No fit.');fig.savefig(BASE/'li_appendix_single22_field_compare.png',dpi=160)
    arrays['radii']=ref['radii'];arrays['phi']=ref['phi'];np.savez_compressed(BASE/'li_appendix_single22_field_compare.npz',**arrays)
    assert source_fingerprint()==fp
    result=dict(status='single22_counterfactual_fixed_Li_raster_comparison',rows=rows,totals=totals,input_sha256=inputs,implementation_sha256=sha(__file__),source_provenance=fp,
      field_normalization='original alpha^-3 / sqrt(cloud mass), unchanged',comparison='9 fixed radii x 72 fixed angles; no amplitude or phase fit; D_RMS uses nominal raster values, intervals retained separately',
      limitations=['Only scalar22 channel was changed; this is not a complete implementation of the printed Appendix for all modes.',
        'No inference that authors used the printed expression in their calculations.',
        'Original finite L6/q12/nr8/h32 and exact-synchronous-vs-a.88 parameter differences retained.',
        'Total flux starts from all unchanged report fluxes and adds the continuous-integral single22 flux difference; raw source quadrature versus continuous interpolation differs slightly.'])
    (BASE/'li_appendix_single22_field_compare.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':main()
