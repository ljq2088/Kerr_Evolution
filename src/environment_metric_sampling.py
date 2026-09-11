"""Parallel radius sampling with resumable, versioned metric tensor caches."""
from concurrent.futures import ProcessPoolExecutor,as_completed
from functools import lru_cache
import hashlib
import importlib.metadata
import json
import multiprocessing
from pathlib import Path
from time import perf_counter
import numpy as np


def initialize(metric,theta):
    global _metric,_theta
    _metric,_theta=metric,theta


def sample_radius(r):
    return np.array([_metric(float(r),float(t)) for t in _theta])


class PrecomputedMetric:
    precomputed=True
    def __init__(self,base,values):
        for name in ('a','r0','m','omega','ellmin','ellmax'):
            setattr(self,name,getattr(base,name))
        self.provenance=base.provenance
        self.values=values

    @lru_cache(maxsize=8192)
    def _values(self,r,theta):
        return tuple(self.values[(r,theta)].ravel())

    def __call__(self,r,theta):
        return np.asarray(self._values(float(r),float(theta))).reshape(4,4).copy()


def precompute_metric(base,radii,theta,cache_root,workers=4):
    if workers<1:raise ValueError('Workers must be positive')
    root=Path(__file__).resolve().parent
    sources=sorted(set(root.glob('lorenz_*.py'))|{root/name for name in (
        'environment_lorenz_mode.py','environment_static_lorenz.py',
        'environment_radial.py','environment_source.py','environment_cloud.py')})
    digest=hashlib.sha256()
    for path in sources:
        digest.update(path.name.encode());digest.update(path.read_bytes())
    metadata=json.dumps(dict(version=1,source_hash=digest.hexdigest(),
        libraries={name:importlib.metadata.version(name) for name in ('pybhpt','numpy','scipy','mpmath')},metric=base.provenance,
        theta=np.asarray(theta).tolist()),sort_keys=True)
    folder=Path(cache_root)/hashlib.sha256(metadata.encode()).hexdigest()[:24]
    folder.mkdir(parents=True,exist_ok=True)
    radii=np.asarray(radii,float);theta=np.asarray(theta,float)
    if len(set(radii))!=len(radii):raise ValueError('Repeated radius')
    values={};missing=[];loaded=0;started=perf_counter()
    def filename(r):return folder/(hashlib.sha256(float(r).hex().encode()).hexdigest()[:24]+'.npz')
    def accept(r,h):
        if h.shape!=(len(theta),4,4) or not np.all(np.isfinite(h)):
            raise ValueError('Invalid cached metric tensor')
        for t,v in zip(theta,h):values[(float(r),float(t))]=v
    for r in radii:
        path=filename(r)
        if path.exists():
            with np.load(path,allow_pickle=False) as saved:
                if str(saved['metadata'])!=metadata or float(saved['r'])!=r:
                    raise ValueError('Metric cache provenance mismatch')
                accept(r,saved['h']);loaded+=1
        else:missing.append(float(r))
    if missing:
        with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn'),
                                 initializer=initialize,initargs=(base,theta)) as executor:
            pending={executor.submit(sample_radius,r):r for r in missing}
            for done,future in enumerate(as_completed(pending),1):
                r=pending[future];h=future.result();accept(r,h)
                path=filename(r);temporary=path.with_suffix('.tmp')
                with temporary.open('wb') as stream:
                    np.savez_compressed(stream,metadata=metadata,r=r,h=h)
                temporary.replace(path)
                print(f'Metric radii {loaded+done}/{len(radii)} saved',flush=True)
    return PrecomputedMetric(base,values),dict(workers=workers,cached_radii=loaded,
        reconstructed_radii=len(missing),elapsed_seconds=perf_counter()-started,
        cache_key=folder.name,source_hash=digest.hexdigest())
