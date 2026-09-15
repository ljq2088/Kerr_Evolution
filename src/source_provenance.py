"""Bind reusable source samples to the implementation that generated them.

Historical reports remain readable. Producers refuse to resume or reuse
unversioned samples rather than relabel old numbers with a new code hash.
"""
import ast,hashlib,importlib.metadata,json,platform
from pathlib import Path


def local_dependency_hashes(root,entrypoints):
    root=Path(root);pending=list(entrypoints);seen={}
    while pending:
        module=pending.pop()
        if module in seen:continue
        path=root/(module+'.py')
        if not path.is_file():continue
        raw=path.read_bytes();seen[module]=hashlib.sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw)):
            if isinstance(node,ast.Import):names=[a.name for a in node.names]
            elif isinstance(node,ast.ImportFrom):names=[node.module or '']
            else:continue
            pending.extend(name.split('.')[0] for name in names)
    return dict(sorted(seen.items()))


def source_fingerprint():
    from lorenz_jet import Jet
    import lorenz_metric
    root=Path(__file__).resolve().parent
    # Include alternate source constructors because run() selects them at runtime.
    files=local_dependency_hashes(root,['report_environment_forced_mode','source_provenance'])
    active=lambda f:f.__module__+'.'+f.__qualname__
    result=dict(schema=1,files=files,python=platform.python_version(),
        libraries={name:importlib.metadata.version(name) for name in ('pybhpt','numpy','scipy','mpmath')},
        runtime=dict(angular=active(lorenz_metric.SpinWeightedSpheroidalHarmonic),
                     kappa=active(lorenz_metric.kappa_jet),jet_dtype=str(Jet.coefficient_dtype)))
    result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True).encode()).hexdigest()
    return result


def samples_hash(samples):
    return hashlib.sha256(json.dumps(samples,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def validate_saved_samples(record,expected):
    if not record.get('samples'):return
    if 'source_provenance' not in record:
        raise ValueError('Unversioned source samples cannot be resumed or reused. '
                         'Keep the historical report and use --output for a fresh calculation.')
    if record['source_provenance']!=expected:
        raise ValueError('Source implementation provenance differs; use --output for a fresh calculation.')
    if record.get('source_samples_sha256')!=samples_hash(record['samples']):
        raise ValueError('Saved source sample checksum differs')
