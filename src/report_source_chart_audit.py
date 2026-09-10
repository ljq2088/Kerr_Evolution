"""Compare scalar contractions in BL and regular advanced Kerr coordinates."""
import json
from pathlib import Path
import numpy as np
from environment_source import ThresholdCloud,kerr_metric
from environment_lorenz_mode import LorenzMetricMode
from environment_ingoing import ingoing_lorenz_source,ingoing_metric


def main():
    cloud=ThresholdCloud(alpha=.3)
    metric=LorenzMetricMode(20.,cloud.a,2,4)
    cases=[]
    for offset in (.05,.005,.0005):
        r=cloud.rp+offset
        h=metric(r,1.1)
        bl=cloud.lorenz_source(r,1.1,h)
        regular=ingoing_lorenz_source(cloud,r,1.1,h)
        case=dict(offset=offset,bl=[float(bl.real),float(bl.imag)],
                  ingoing=[float(regular.real),float(regular.imag)],
                  relative_difference=float(abs(bl-regular)/abs(regular)),
                  bl_metric_condition=float(np.linalg.cond(kerr_metric(r,1.1,cloud.a))),
                  ingoing_metric_condition=float(np.linalg.cond(ingoing_metric(r,1.1,cloud.a))))
        cases.append(case)
        print(case,flush=True)
    out=Path(__file__).resolve().parents[1]/'docs/environment_reproduction/source_chart_audit.json'
    out.write_text(json.dumps(dict(status='coordinate_cross_check_not_flux_convergence',cases=cases),indent=2)+'\n')


if __name__=='__main__':
    main()
