"""User-scoped continuation: finish Fig.2 at r0=18.3 only, then exit.
Never launch the full Fig.2 scan or Fig.9/Fig.10 follow-up.
"""
import os,sys,json,subprocess,traceback,importlib.util
from pathlib import Path
OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
spec=importlib.util.spec_from_file_location("fig2_scoped_backend",OUT/"run_fig2.py")
b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
b.ORBIT_PLAN=[20.,18.3]
def main():
    fp=b.fingerprint()
    old=json.loads((OUT/"execution.json").read_text())
    if old["implementation_sha256"]!=fp:raise RuntimeError("Implementation changed; refuse mixed continuation")
    state=dict(status="running",pid=os.getpid(),started_utc=b.stamp(),implementation_sha256=fp,orbit_plan=[18.3],workers=4,completed_groups=[],failures=[],scope="Finish only Fig.2 r0=18.3M then stop; full scan and Fig.9/10 continuation cancelled by user")
    seq=dict(status="running_fig2_r18p3_only",pid=os.getpid(),order=["Fig2 r0=18.3M only"],workers=4,auto_continue_fig9_fig10=False)
    b.save(OUT/"sequence.json",seq);b.save(OUT/"execution.json",state)
    try:
        b.ingest()
        for mg in [6,1,2,3,4,5]:
            if all(b.modepath(18.3,*t).exists() for t in b.targets(18.3,mg)):continue
            if b.fingerprint()!=fp:raise RuntimeError("Numerical implementation changed")
            state.update(active_orbit=18.3,active_metric_abs_m=mg,updated_utc=b.stamp());b.save(OUT/"execution.json",state)
            with (OUT/f"group_r18.3_mg{mg}.log").open("a") as log:
                code=subprocess.call([sys.executable,"-u",str(OUT/"run_fig2.py"),"--group","18.3",str(mg)],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
            if code:raise RuntimeError(f"Group {mg} exited {code}; preserved checkpoints")
            state["completed_groups"].append([18.3,mg]);b.save(OUT/"execution.json",state);b.render()
        missing=[(l,m) for l,m in b.EXPECTED if not b.modepath(18.3,1,l,m).exists()]
        if missing:raise RuntimeError(f"Incomplete requested radius: {missing}")
        for l,m in b.EXPECTED:b.check_hashes(json.loads(b.modepath(18.3,1,l,m).read_text())["implementation_sha256"])
        b.render();comparison=json.loads((OUT/"comparison.json").read_text())
        summary={"status":"completed_requested_radius_finite_resolution_convergence_pending","r0":18.3,"completed_modes":len(b.EXPECTED),"comparison":[x for x in comparison["rows"] if x["r0"]==18.3],"completed_utc":b.stamp(),"scope":"Only r0=18.3M requested; no further figure jobs launched"}
        b.save(OUT/"r18p3_completion.json",summary)
        state.update(status=summary["status"],completed_utc=b.stamp(),pid=None)
        seq.update(status="completed_requested_radius_stopped",completed_utc=b.stamp(),pid=None)
    except Exception:
        state.update(status="failed_resumable",traceback=traceback.format_exc(),pid=None)
        seq.update(status="needs_attention",pid=None)
        raise
    finally:
        b.save(OUT/"execution.json",state);b.save(OUT/"sequence.json",seq)
if __name__=="__main__":main()
