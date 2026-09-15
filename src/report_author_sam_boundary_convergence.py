"""Compare two separately completed original-author boundary-control runs."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from report_author_hdf5_alignment import dec,errors

def main():
    p=argparse.ArgumentParser();p.add_argument('--coarse',required=True);p.add_argument('--fine',required=True);p.add_argument('--output',required=True);args=p.parse_args()
    paths=[Path(args.coarse),Path(args.fine)];coarse,fine=[json.loads(x.read_text()) for x in paths]
    for k in ['status','a','rp','m','lmax','protected_output_ells']:
        if coarse[k]!=fine[k]:raise ValueError('comparison parameter mismatch '+k)
    guard=fine['protected_output_ells'];rows=[]
    for c,f in zip(coarse['rows'],fine['rows'],strict=True):
        if c['r']!=f['r']:raise ValueError('radius mismatch')
        rows.append(dict(r=f['r'],side=f['side'],reference_boundary_change=errors(dec(c['reference'])[guard],dec(f['reference'])[guard]),local_vs_fine=f['errors_protected']))
    result=dict(status='completed_original_author_boundary_refinement',scope='Three-radii original Sam L10 m1, protected spherical ell1..3; not endpoint or high-L convergence and not final flux comparison',a=fine['a'],rp=fine['rp'],m=fine['m'],lmax=fine['lmax'],protected_output_ells=guard,coarse_controls=dict(inford=6,horord=5,rinf=10000,xhor=.0001,precision=50,accuracy_goal=24),fine_controls=dict(inford=7,horord=6,rinf=20000,xhor=.0001,precision=50,accuracy_goal=24),inputs=[dict(path=str(x),sha256=hashlib.sha256(x.read_bytes()).hexdigest()) for x in paths],rows=rows)
    result['maximum_relative_author_boundary_change']=max(e['relative_l2'] for r in rows for e in r['reference_boundary_change'] if e['relative_l2'] is not None)
    result['maximum_relative_local_vs_fine']=max(e['relative_l2'] for r in rows for e in r['local_vs_fine'] if e['relative_l2'] is not None)
    Path(args.output).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for row in rows:print(row['r'],max(e['relative_l2'] for e in row['reference_boundary_change'] if e['relative_l2'] is not None),max(e['relative_l2'] for e in row['local_vs_fine'] if e['relative_l2'] is not None))
if __name__=='__main__':main()
