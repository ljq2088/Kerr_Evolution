"""Freeze a verifiable Git snapshot without stopping or rewriting a live run.

Only modes listed in one atomic execution.json read are selected. Mode, radial
and metric-bank files are immutable after completion. Logs are copied at their
observed length; the running process continues to append to its own log.
"""
import argparse,hashlib,json
from pathlib import Path
from datetime import datetime,timezone

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--run',type=Path,default=Path('docs/figure1_alignment_20260926'))
    ap.add_argument('--tag',required=True)
    args=ap.parse_args()
    if not args.tag.replace('_','').replace('-','').isalnum():raise ValueError('Invalid snapshot tag')
    folder=args.run
    state=json.loads((folder/'execution.json').read_text())
    target=folder/f'snapshot_{args.tag}.json'
    if target.exists():raise FileExistsError('Snapshots must not be overwritten')
    selected={}
    modes=state['completed_modes']
    if len(set(map(tuple,modes)))!=len(modes):raise ValueError('Duplicate completed mode')
    for ell,m in modes:
        p=folder/f'mode_l{ell}_m{m}.json';d=json.loads(p.read_text())
        if d['config']!=state['config'] or d['implementation_sha256']!=state['implementation_sha256']:
            raise ValueError('Mixed mode provenance')
        radial=folder/d['radial_file'];bank=folder/f'metric_mg{m-1}.npz'
        if digest(radial)!=d['radial_sha256'] or digest(bank)!=d['metric_bank_sha256']:
            raise ValueError('Mode dependency checksum mismatch')
        for file in (p,radial,bank):selected[file.name]=digest(file)
    if (folder/'static_matching.json').exists():
        selected['static_matching.json']=digest(folder/'static_matching.json')
    for log in ('run.log','static_matching.log'):
        source=folder/log
        if source.exists():
            output=folder/f'{source.stem}_{args.tag}.log'
            if output.exists():raise FileExistsError(output)
            output.write_bytes(source.read_bytes())
            selected[output.name]=digest(output)
    result=dict(snapshot_utc=datetime.now(timezone.utc).isoformat(),snapshot_kind='partial_immutable_publication',
        completed_count=len(modes),expected_count=state['expected_mode_count'],
        completed_modes=modes,execution_at_snapshot=state,files_sha256=selected,
        caution='Remote snapshot is not a live process and is not a full field result. Missing modes remain missing.')
    target.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(dict(snapshot=str(target),completed=len(modes),files=len(selected))))
if __name__=='__main__':main()
