"""Sensitivity to decimal precision of the published cloud frequency.
Sources remain frozen: this is NOT a new physical background or a fit.
"""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.optimize import minimize_scalar
from environment_radial import RadialGreen
from environment_response import SampledResponse
from report_li_fixed_source_background import OriginalSourceCoordinate
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'docs/root_cause_followup_20260917'
def pair(z):return [float(z.real),float(z.imag)]
def main():
    rows=[];inputs={}
    for orbit in [41.1,42.1]:
        f=ROOT/f'docs/field_alignment_20260917/nonstatic/rp{str(orbit).replace(".","p")}_sl2_sm2_L6_q12_nr8_h32.json'
        data=json.loads(f.read_text());p=data['parameters'];inputs[str(f.relative_to(ROOT))]=hashlib.sha256(f.read_bytes()).hexdigest()
        nodes=np.array([v['r'] for v in data['samples']]);J=np.array([complex(*v['source']) for v in data['samples']])
        oldrp=1+np.sqrt(1-p['metric']['a']**2)
        rr=np.array([20.,50.,70.,100.,120.,150.,180.]);grid=np.linspace(5,200,1200)
        baseline=None
        for label,wc in [('exact_Re',.2962935347256114628),('printed_6_decimals',.296294),('rounded_5_decimals',.29629),('lower_last_printed_digit',.2962935),('upper_last_printed_digit',.2962945)]:
            a=.88;w=wc+1/(orbit**1.5+a)
            g=RadialGreen(a,.3,w,2,2,rmax=32000,offset=1e-4,rtol=1e-11,infinity_method='coulomb')
            response=SampledResponse(OriginalSourceCoordinate(g,oldrp),p['source_panels'],nodes,J,log_first=True)
            value=response.evaluate(rr)[0];vals=abs(response.evaluate(grid)[0]);mins=[]
            for i in range(1,len(grid)-1):
                if vals[i]<vals[i-1] and vals[i]<vals[i+1]:
                    x=minimize_scalar(lambda x:float(abs(response.evaluate([x])[0][0])**2),bounds=(grid[i-1],grid[i+1]),method='bounded',options={'xatol':1e-9})
                    mins.append(dict(r=float(x.x),amplitude=float(np.sqrt(x.fun))))
            if baseline is None:baseline=value.copy();zh0=response.horizon_coefficient
            rows.append(dict(orbit=orbit,label=label,omega_cloud=wc,forced_omega=w,radii=rr.tolist(),field=[pair(z) for z in value],
                ratio_to_exact_Re=[pair(z) for z in value/baseline],max_relative_field_change=float(np.max(abs(value/baseline-1))),
                sampled_field_relative_L2=float(np.linalg.norm(value-baseline)/np.linalg.norm(baseline)),
                relative_ZH_change=pair(response.horizon_coefficient/zh0-1),radial_minima=mins))
            print(orbit,label,rows[-1]['max_relative_field_change'],mins,flush=True)
    result=dict(status='fixed_source_decimal_frequency_sensitivity_not_consistent_new_background',
        reason='Li Eq omega_c_11 displays 0.296294 without publishing the full-precision frequency used for its forced-field plot. No parameter is optimized.',
        limitations=['Only radial Green frequency and spin change. Cloud, metric and projected source remain frozen at synchronized baseline.',
        'Printed decimal rounding is not proof the author truncated the actual computation.'],rows=rows,input_sha256=inputs,
        implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/'scalar22_frequency_rounding_control.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
