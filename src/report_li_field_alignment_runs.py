"""Fresh finite-resolution Li Fig.11 nonstatic channels, with metric reuse.

This diagnostic orchestrator changes no production equations or normalizations.
Each orbit is a separate process with four radius workers; positive metric m is
sampled once, and only its real-metric Fourier conjugate is reused for negative m.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "docs/field_alignment_20260917"
OUTPUT = BASE / "nonstatic"
POSITIVE_SIX = [(2,2), (3,3), (4,2), (4,4), (5,3), (5,5)]

def stamp():
    return datetime.now(timezone.utc).isoformat()

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)

def orbit_tag(rp):
    return f"{rp:g}".replace(".", "p")

def channel_path(rp, ell, m):
    return OUTPUT / f"rp{orbit_tag(rp)}_sl{ell}_sm{m}_L6_q12_nr8_h32.json"

def channels(mg):
    return [(ell,m) for ell in range(2,6) for m in range(-ell,ell+1,2)
            if m != 1 and abs(m-1) == mg]

def make_args(rp, ell, m, method=None):
    from report_environment_forced_mode import argument_parser
    near = m == 2
    method = method or ("coulomb" if near else "series")
    return argument_parser().parse_args([
        "--output", str(channel_path(rp,ell,m)),
        "--alpha", ".3", "--orbital-radius", str(rp),
        "--metric-ellmax", "6", "--metric-m", str(m-1),
        "--scalar-ell", str(ell), "--angular-order", "12",
        "--radial-order", "8", "--horizon-order", "32", "--horizon-log",
        "--source-inner-offset", ".0005", "--source-outer-radius", "320",
        "--green-outer-radius", "32000" if near else "4000",
        "--infinity-method", method])

def worker(rp, workers):
    from environment_source import ThresholdCloud
    from environment_lorenz_mode import LorenzMetricMode, ConjugateMetricMode
    from environment_metric_sampling import precompute_metric
    from environment_dense_metric import configure_metric_backend
    from environment_radial import RadialGreen
    from report_environment_forced_mode import source_grid, run
    from source_provenance import source_fingerprint, validate_saved_samples
    configure_metric_backend(False)
    cloud = ThresholdCloud(alpha=.3)
    fingerprint = source_fingerprint()
    path = BASE / f"nonstatic_rp{orbit_tag(rp)}_manifest.json"
    record = dict(status="running",started_utc=stamp(),pid=os.getpid(),
        implementation_sha256=sha(__file__),source_provenance=fingerprint,
        parameters=dict(alpha=.3,orbital_radius=rp,a=cloud.a,
          cloud_omega=cloud.omega,stationary_synchronized_cloud=True,
          Li_label_a=.88,metric_ellmax=6,angular_order=12,radial_order=8,
          horizon_order=32,source_inner_offset=.0005,source_outer_radius=320,
          workers=workers,positive_six=[list(v) for v in POSITIVE_SIX]),
        groups=[],channels=[],failures=[],
        limitations=["Finite L6/q12/nr8/h32 result; no convergence claim.",
          "Stationary cloud spin differs from Li figure label a=.88.",
          "No static metric mg0 channels are calculated here.",
          "Both signs are retained as separately projected scalar channels; only the real metric Fourier coefficients are conjugated.",
          "Inverse-r series diagnostic is a necessary local boundary check, not a full cutoff-convergence proof."])
    save(path,record)
    try:
        for mg in range(1,7):
            cases=sorted(channels(mg),key=lambda x:(x not in POSITIVE_SIX, x))
            if not cases:
                continue
            settings={}
            for ell,m in cases:
                args=make_args(rp,ell,m)
                omega=cloud.omega+(m-1)/(rp**1.5+cloud.a)
                started=time.perf_counter()
                green=RadialGreen(cloud.a,cloud.mu,omega,ell,m,
                    rmax=args.green_outer_radius,offset=args.green_horizon_offset,
                    rtol=1e-11,infinity_method=args.infinity_method)
                alternate=None
                if args.infinity_method=="series" and green.series_last_term_relative>1e-3:
                    alternate=dict(rejected_method="series",
                        series_last_term_relative=green.series_last_term_relative,
                        reason="Production unresolved-series guard")
                    args=make_args(rp,ell,m,"coulomb")
                    green=RadialGreen(cloud.a,cloud.mu,omega,ell,m,
                        rmax=args.green_outer_radius,offset=args.green_horizon_offset,
                        rtol=1e-11,infinity_method=args.infinity_method)
                probes=np.array([cloud.rp+.0005,3.,20.,rp,100.,320.])
                spread=float(np.max(abs(green.wronskian(probes)/green.w0-1)))
                if not np.isfinite(spread) or spread>1e-6:
                    raise ValueError(f"Unacceptable Green Wronskian spread {spread} for {(ell,m)}")
                settings[(ell,m)]=(args,dict(method=green.infinity_method,
                    outer_radius=green.rmax,propagating=bool(green.propagating),
                    series_last_term_relative=green.series_last_term_relative,
                    wronskian_relative_spread=spread,alternate=alternate,
                    elapsed_seconds=time.perf_counter()-started))
            base=LorenzMetricMode(rp,cloud.a,mg,6)
            gridargs=next(iter(settings.values()))[0]
            _,radii,_=source_grid(gridargs,cloud)
            theta=np.arccos(np.polynomial.legendre.leggauss(12)[0])
            record["active_metric_m"]=mg
            record["active_metric_started_utc"]=stamp()
            save(path,record)
            print(f"ORBIT {rp:g} MG {mg} PRECOMPUTE START ({len(radii)} radii)",flush=True)
            sampled,audit=precompute_metric(base,radii,theta,ROOT/"outputs/metric_cache",workers)
            if source_fingerprint()!=fingerprint:
                raise RuntimeError("Production source changed during metric sampling")
            record["groups"].append(dict(metric_m=mg,precomputation=audit,
                completed_utc=stamp(),source_provenance=fingerprint))
            save(path,record)
            for ell,m in cases:
                args,boundary=settings[(ell,m)]
                chosen=sampled if m-1>0 else ConjugateMetricMode(sampled)
                print(f"ORBIT {rp:g} SCALAR {(ell,m)} START",flush=True)
                before=time.perf_counter()
                out,result=run(args,cloud=cloud,metric=chosen)
                validate_saved_samples(result,fingerprint)
                if source_fingerprint()!=fingerprint:
                    raise RuntimeError("Production source changed during scalar projection")
                result["field_alignment_audit"]=dict(
                    run_manifest=str(path.relative_to(ROOT)),generated_utc=stamp(),
                    positive_six=(ell,m) in POSITIVE_SIX,
                    metric_sampling=audit,metric_positive_m=mg,
                    metric_negative_m_obtained_by_conjugation=m-1<0,
                    boundary_preflight=boundary,
                    finite_resolution=True,Li_label_a=.88,
                    exact_a_difference=cloud.a-.88,
                    no_fitted_amplitude_or_phase=True)
                save(out,result)
                row=dict(scalar_ell=ell,scalar_m=m,metric_m=m-1,
                    positive_six=(ell,m) in POSITIVE_SIX,
                    file=str(out.relative_to(ROOT)),sha256=sha(out),
                    source_samples_sha256=result["source_samples_sha256"],
                    z_h=result["z_h"],z_inf=result["z_inf"],
                    flux=result["flux"],
                    elapsed_seconds=time.perf_counter()-before,
                    wronskian_relative_spread=result["wronskian_relative_spread"])
                record["channels"].append(row)
                save(path,record)
                print(f"ORBIT {rp:g} SCALAR {(ell,m)} COMPLETE",flush=True)
        record["status"]="finite_resolution_nonstatic_channels_completed"
        record["completed_utc"]=stamp()
        record.pop("active_metric_m",None)
        save(path,record)
        return 0
    except Exception as exc:
        record["status"]="failed_resumable"
        record["failures"].append(dict(time_utc=stamp(),exception=repr(exc),
            traceback=traceback.format_exc()))
        save(path,record)
        raise

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--worker-rp",type=float)
    parser.add_argument("--workers",type=int,default=4)
    args=parser.parse_args()
    BASE.mkdir(parents=True,exist_ok=True);OUTPUT.mkdir(parents=True,exist_ok=True)
    if args.worker_rp is not None:
        return worker(args.worker_rp,args.workers)
    record=dict(status="running",started_utc=stamp(),orbits=[41.1,42.1],
        metric_groups_parallel=2,metric_radius_workers_per_group=args.workers,
        output_directory=str(OUTPUT.relative_to(ROOT)),results=[],
        implementation_sha256=sha(__file__))
    manifest=BASE/"nonstatic_execution_manifest.json"
    save(manifest,record)
    def launch(rp):
        log=BASE/f"nonstatic_rp{orbit_tag(rp)}.log"
        command=[sys.executable,"-u",str(Path(__file__).resolve()),
                 "--worker-rp",str(rp),"--workers",str(args.workers)]
        with log.open("a",buffering=1) as stream:
            stream.write(f"\nLAUNCH {stamp()}\n")
            process=subprocess.run(command,cwd=ROOT,stdout=stream,
                stderr=subprocess.STDOUT,env=dict(os.environ,OPENBLAS_NUM_THREADS="1",
                OMP_NUM_THREADS="1",MKL_NUM_THREADS="1",PYTHONPATH=str(ROOT/"src")))
        return dict(orbital_radius=rp,exit_code=process.returncode,
            log=str(log.relative_to(ROOT)),completed_utc=stamp(),command=command)
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures=[pool.submit(launch,rp) for rp in record["orbits"]]
        for future in as_completed(futures):
            row=future.result();record["results"].append(row);save(manifest,record)
            print(json.dumps(row),flush=True)
    record["status"]=("finite_resolution_nonstatic_channels_completed"
        if all(v["exit_code"]==0 for v in record["results"]) else "failed_resumable")
    record["completed_utc"]=stamp();save(manifest,record)
    return int(record["status"]=="failed_resumable")

if __name__=="__main__":
    raise SystemExit(main())
