"""Check the scalar source is invariant under independent complex radial rescaling."""
import hashlib, importlib, json, time
from pathlib import Path
import numpy as np
from pybhpt.radial import RadialTeukolsky as NativeRadial
from source_provenance import local_dependency_hashes
ROOT=Path(__file__).resolve().parents[1]

def main():
    baseline_path=ROOT/'docs/environment_reproduction/metric_backend_response_L1_s2AUTO_gAUTO_rtolNone.json'
    baseline=json.loads(baseline_path.read_text())
    for name,expected in baseline['implementation_sha256'].items():
        assert hashlib.sha256((ROOT/'src'/f'{name}.py').read_bytes()).hexdigest()==expected,name
    for item in baseline['package_artifacts'].values():
        assert hashlib.sha256(Path(item['path']).read_bytes()).hexdigest()==item['sha256']
    import lorenz_metric
    class ScaledRadial(NativeRadial):
        def scale(self,bc):
            if self.spinweight not in (0,-1,1): return 1.
            side=1 if bc=='Up' else 0
            return (1+.3*(self.spinweight+2)+.4*side)*np.exp(1j*(.4*self.spinweight+.3*self.azimuthalmode+.8*side))
        def radialsolution(self,bc,i):
            return self.scale(bc)*super().radialsolution(bc,i)
        def radialderivative(self,bc,i):
            return self.scale(bc)*super().radialderivative(bc,i)
        def radialsolutions(self,bc):
            return self.scale(bc)*super().radialsolutions(bc)
        def radialderivatives(self,bc):
            return self.scale(bc)*super().radialderivatives(bc)
    for name in ('lorenz_metric','lorenz_spin1','lorenz_spin1_chiral','lorenz_chi'):
        module=importlib.import_module(name)
        for obj in vars(module).values():
            if callable(getattr(obj,'cache_clear',None)):obj.cache_clear()
        module.RadialTeukolsky=ScaledRadial
    from environment_source import ThresholdCloud,project_source
    from environment_lorenz_mode import LorenzMetricMode,ConjugateMetricMode
    from environment_angular_diagnostic import install_dense_angular_diagnostic
    install_dense_angular_diagnostic()
    cloud=ThresholdCloud();metric=ConjugateMetricMode(LorenzMetricMode(20.,cloud.a,1,1))
    radii=np.array([s['r'] for s in baseline['samples']])
    indices=sorted(set(int(np.argmin(abs(radii-r))) for r in [1.481,2.,10.,19.9,20.1,40.,120.,316.]))
    indices=sorted(set(indices+[int(np.flatnonzero(radii>20)[0])]))
    output=ROOT/'docs/environment_reproduction/gauge_radial_normalization_audit.json'
    result=dict(status='sampling',scaling='For s=0,+/-1, C_In/Up=(1+.3*(s+2)+.4*side)*exp(i*(.4*s+.3*m+.8*side))',
                baseline_sha256=hashlib.sha256(baseline_path.read_bytes()).hexdigest(),
                implementation_sha256=local_dependency_hashes(ROOT/'src',['report_gauge_radial_normalization_audit']),
                scope='Metric ell1 -> scalar00 source, all s0 and +/-1 branches including negative chirality',rows=[])
    started=time.perf_counter()
    for index in indices:
        r=radii[index];_,new=project_source(cloud,[r],20.,0,0,metric,ntheta=18)
        old=complex(*baseline['samples'][index]['source']);value=new[0]
        result['rows'].append(dict(r=float(r),source_original=[old.real,old.imag],source_rescaled=[value.real,value.imag],relative_difference=float(abs(value-old)/abs(old))))
        print(r,result['rows'][-1]['relative_difference'],flush=True)
    result.update(status='completed',max_source_relative_difference=max(x['relative_difference'] for x in result['rows']),elapsed_seconds=time.perf_counter()-started,
                  limitations=['Normalization invariance checks consistent use of radial amplitudes; it cannot establish correct physical source coefficients','Finite eight-node ell1 diagnostic, not a new full flux calculation'])
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':main()
