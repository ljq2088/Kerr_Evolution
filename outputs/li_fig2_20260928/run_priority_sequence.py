import os,sys,json,subprocess,shutil,traceback
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent;NEXT=ROOT/'outputs/li_fig9_10_20260928'
def save(d):
 d['updated_utc']=datetime.now(timezone.utc).isoformat();tmp=OUT/'sequence.tmp';tmp.write_text(json.dumps(d,indent=2));tmp.replace(OUT/'sequence.json')
state={'status':'running_fig2','pid':os.getpid(),'order':['Fig2','Fig10'],'workers':4};save(state)
try:
 with (OUT/'run_priority.log').open('a') as log:
  code=subprocess.call([sys.executable,'-u',str(OUT/'run_fig2.py'),'--run'],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
 d=json.loads((OUT/'execution.json').read_text())
 if code or d['status']!='completed_finite_scan_convergence_pending':
  state.update(status='fig2_needs_attention',returncode=code);save(state);sys.exit(1)
 state['status']='reusing_fig2_results';save(state)
 # Both jobs are stopped at this transition; publish only committed cache files.
 cache=ROOT/'outputs/metric_cache';copied=0
 for f in (OUT/'metric_cache').glob('*/*.npz'):
  dest=cache/f.parent.name/f.name
  if not dest.exists():
   dest.parent.mkdir(parents=True,exist_ok=True);tmp=dest.with_suffix('.transfer');shutil.copy2(f,tmp);tmp.replace(dest);copied+=1
 for f in OUT.glob('r*_c1_l*_m*.json'):
  dest=NEXT/f.name
  if not dest.exists():
   tmp=dest.with_suffix('.transfer');shutil.copy2(f,tmp);tmp.replace(dest)
 state.update(status='running_fig10',metric_cache_files_reused=copied);save(state)
 with (NEXT/'run_after_fig2.log').open('a') as log:
  code=subprocess.call([sys.executable,'-u',str(NEXT/'run_comparison.py'),'--run','--workers','4'],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
 final=json.loads((NEXT/'execution.json').read_text())
 state.update(status='completed_finite_scans_convergence_pending' if code==0 and final['status']=='completed_finite_scan_convergence_pending' else 'fig10_needs_attention',returncode=code);save(state)
except Exception:
 state.update(status='failed_resumable',traceback=traceback.format_exc());save(state);raise
