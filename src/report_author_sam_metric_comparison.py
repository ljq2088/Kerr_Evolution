"""Unfitted original Sam Dolan run vs local metric sector comparison."""
import argparse,json,hashlib,re
from pathlib import Path
import numpy as np
from report_author_hdf5_alignment import enc,dec,errors
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read_wolfram_json(path):
    # Wolfram 14 RawJSON exports exact numerical zeros with uncertainty as
    # 0.e-49, which is not strict JSON. Repair decimal syntax outside strings
    # only; preserve the original author file and its SHA256 unchanged.
    repairs=[]
    def replace(match):
        if match.group(1) is not None:return match.group(0)
        repairs.append(match.group(0))
        return match.group(2)+'.0'+match.group(3)
    pattern=r'("(?:\\.|[^"\\])*")|(-?\d+)\.([eE][+-]?\d+)'
    raw=Path(path).read_text();normalized=re.sub(pattern,replace,raw)
    return json.loads(normalized),dict(decimal_exponent_repairs=len(repairs),unique_tokens=sorted(set(repairs)),policy='Outside-string decimal syntax only; original bytes and numeric values retained')

def main():
    p=argparse.ArgumentParser();p.add_argument('--reference',required=True);p.add_argument('--local',nargs='+',required=True);p.add_argument('--output',required=True);args=p.parse_args()
    rp=Path(args.reference);lps=[Path(x) for x in args.local];ref,normalization=read_wolfram_json(rp);locals_=[json.loads(lp.read_text()) for lp in lps];local=dict(locals_[0]);local['rows']=sorted([row for d in locals_ for row in d['rows']],key=lambda row:row['r'])
    for d in locals_:
        if d['status']!='completed' or any(d[k]!=local[k] for k in ['a','rp','m','lmax','q']):raise ValueError('local sampler mismatch/incomplete')
    if local['status']!='completed':raise ValueError('local run incomplete')
    for k in ['a','rp','m','lmax']:
        if abs(ref[k]-local[k])>1e-13:raise ValueError('parameter mismatch '+k)
    L=local['lmax'];protected=np.arange(abs(local['m']),L-6)
    result=dict(status='completed_original_author_run_comparison',scope='Original Sam Dolan radiative code, independently executed at matched parameters; not original 2025 environmental-flux production data',a=local['a'],rp=local['rp'],m=local['m'],lmax=L,local_q=local['q'],protected_output_ells=protected.tolist(),protection='Strict L-ell>6 guard from paper; finite cutoff still needs convergence',reference_path=str(rp),reference_sha256=sha(rp),local_inputs=[dict(path=str(lp),sha256=sha(lp)) for lp in lps],reference_json_normalization=normalization,basis='All compared components spin-weighted spherical, trace explicitly transformed by original author Bmat0',scale_phase='No fitting or substitution',rows=[])
    for rr,lr in zip(ref['rows'],local['rows'],strict=True):
        if not rr['all_numeric'] or abs(rr['r']-lr['r'])>1e-13:raise ValueError('incomplete/mismatched sample')
        rs=np.array([dec(rr[k]).T for k in ['spin0','spin1','spin2']]);rs[0,:,9]=dec(rr['trace_spherical'])
        ls=dec(lr['local_sectors']);rt=rs.sum(axis=0);lt=ls.sum(axis=0)
        result['rows'].append(dict(r=rr['r'],side=rr['side'],reference_sectors=enc(rs),local_sectors=enc(ls),reference=enc(rt),local=enc(lt),errors_all_output_ells=errors(lt,rt),errors_protected=errors(lt[protected],rt[protected]),errors_protected_by_sector=[dict(sector=i,components=errors(ls[i,protected],rs[i,protected])) for i in range(3)]))
    Path(args.output).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for row in result['rows']:
        print(row['side'],'total protected errors',[e['relative_l2'] for e in row['errors_protected']])
        for sector in row['errors_protected_by_sector']:
            print('sector',sector['sector'],[e['relative_l2'] for e in sector['components']])
if __name__=='__main__':main()
