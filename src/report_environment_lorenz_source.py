"""First actual Kerr/cloud source projection, explicitly not a flux reproduction."""
import json
from pathlib import Path
from environment_source import ThresholdCloud,project_source
from environment_lorenz_mode import LorenzMetricMode


def main():
    cloud=ThresholdCloud(alpha=.3)
    metric=LorenzMetricMode(20.,cloud.a,2,4)
    cases=[]
    for ntheta in (6,10):
        omega,source=project_source(cloud,[10.],20.,3,3,metric,ntheta=ntheta)
        case=dict(ntheta=ntheta,omega=float(omega),source=[float(source[0].real),float(source[0].imag)])
        cases.append(case)
        print(case,flush=True)
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/actual_lorenz_source_development.json'
    out.write_text(json.dumps(dict(status='local_source_only_not_a_waveform_or_flux',metric=metric.provenance,
                   alpha=.3,cloud_mass=1.,r=10.,scalar_ell=3,scalar_m=3,cases=cases),indent=2)+'\n')


if __name__=='__main__':
    main()
