import copy
from pathlib import Path
import sys
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from source_provenance import local_dependency_hashes,samples_hash,validate_saved_samples


def test_dependency_fingerprint_covers_transitive_local_imports(tmp_path):
    (tmp_path/'driver.py').write_text('from physics import answer\nimport numpy\n')
    (tmp_path/'physics.py').write_text('from coefficients import value\n')
    (tmp_path/'coefficients.py').write_text('value=1\n')
    old=local_dependency_hashes(tmp_path,['driver'])
    (tmp_path/'coefficients.py').write_text('value=2\n')
    new=local_dependency_hashes(tmp_path,['driver'])
    assert set(old)=={'driver','physics','coefficients'}
    assert old['driver']==new['driver'] and old['coefficients']!=new['coefficients']


def test_reuse_rejects_legacy_changed_code_and_modified_samples():
    rows=[dict(r=3.,weight=.2,source=[1.,2.])];provenance={'schema':1,'sha256':'current'}
    legacy=dict(samples=rows)
    with pytest.raises(ValueError,match='Unversioned'):validate_saved_samples(legacy,provenance)
    valid=dict(samples=rows,source_provenance=provenance,source_samples_sha256=samples_hash(rows))
    validate_saved_samples(valid,provenance)
    with pytest.raises(ValueError,match='provenance differs'):validate_saved_samples(valid,{'schema':1,'sha256':'changed'})
    corrupt=copy.deepcopy(valid);corrupt['samples'][0]['source'][0]=3.
    with pytest.raises(ValueError,match='checksum'):validate_saved_samples(corrupt,provenance)
    validate_saved_samples({'samples':[]},provenance)


def test_forced_run_resume_uses_version_guard_before_reusing_samples(tmp_path,monkeypatch):
    import json
    import numpy as np
    import report_environment_forced_mode as forced
    from environment_source import ThresholdCloud
    # A deterministic test source exercises persistence without reconstructing
    # a costly metric; this is an I/O regression, not a physics validation.
    cloud=ThresholdCloud();metric=type('TestMetric',(),{})()
    metric.r0=20.;metric.a=cloud.a;metric.m=-1;metric.omega=-1/(20.**1.5+cloud.a)
    metric.ellmax=metric.ellmin=1
    metric.provenance=dict(a=cloud.a,m_g=-1,ellmax=1,orbital_radius=20.)
    monkeypatch.setattr(forced,'source_fingerprint',lambda:{'schema':1,'sha256':'test-current'})
    calls=[]
    def source(cloud,radii,r0,ell,m,metric,ntheta):
        calls.extend(radii)
        return cloud.omega+metric.omega,np.exp(-np.asarray(radii)/10)*1e-5
    monkeypatch.setattr(forced,'project_source',source)
    output=tmp_path/'response.json'
    args=forced.argument_parser().parse_args(['--metric-m','-1','--scalar-ell','0',
         '--metric-ellmax','1','--angular-order','2','--source-outer-radius','80','--output',str(output)])
    forced.run(args,cloud=cloud,metric=metric)
    count=len(calls);assert count>0
    forced.run(args,cloud=cloud,metric=metric);assert len(calls)==count
    old=json.loads(output.read_text());old.pop('source_provenance');output.write_text(json.dumps(old))
    with pytest.raises(ValueError,match='Unversioned'):
        forced.run(args,cloud=cloud,metric=metric)
    assert len(calls)==count
