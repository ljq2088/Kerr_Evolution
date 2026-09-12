import json
import sys
from pathlib import Path
from types import SimpleNamespace
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
import report_environment_mode_batch as batch


def test_channel_failure_preserves_completed_results_and_records_terminal_state(tmp_path, monkeypatch):
    folder = tmp_path / 'docs' / 'environment_reproduction'
    folder.mkdir(parents=True)
    monkeypatch.setattr(batch, '__file__', str(tmp_path / 'src' / 'report.py'))
    monkeypatch.setattr(sys, 'argv', ['report.py', '--scalar-ells', '3', '5', '--metric-m', '0'])
    cloud = SimpleNamespace(m=1)
    metric = SimpleNamespace(m=0, _values=SimpleNamespace(
        cache_info=lambda: SimpleNamespace(hits=0, misses=0)))
    monkeypatch.setattr(batch, 'build_cloud', lambda args: cloud)
    monkeypatch.setattr(batch, 'build_metric', lambda args, shared: metric)
    calls = []

    def solve(args, **kwargs):
        calls.append(args.scalar_ell)
        if args.scalar_ell == 5:
            raise ValueError('Unresolved infinity boundary')
        return folder / 'first.json', {'flux': {'horizon': {'orbital_energy': 0.}}}

    monkeypatch.setattr(batch, 'run', solve)
    with pytest.raises(ValueError, match='Unresolved infinity boundary'):
        batch.main()
    paths = list(folder.glob('scalar_batch_*.json'))
    assert len(paths) == 1 and calls == [3, 5]
    result = json.loads(paths[0].read_text())
    assert result['status'] == 'batch_failed_partial_results_preserved'
    assert [c['scalar_ell'] for c in result['channels']] == [3]
    assert result['channels'][0]['file'] == 'first.json'
    assert result['failure'] == dict(scalar_ell=5, scalar_m=1,
        error_type='ValueError', message='Unresolved infinity boundary')


def test_metric_precomputation_failure_records_terminal_state(tmp_path, monkeypatch):
    import environment_metric_sampling
    folder = tmp_path / 'docs' / 'environment_reproduction'
    folder.mkdir(parents=True)
    monkeypatch.setattr(batch, '__file__', str(tmp_path / 'src' / 'report.py'))
    monkeypatch.setattr(sys, 'argv', ['report.py', '--scalar-ells', '2', '--metric-m', '1', '--workers', '2'])
    monkeypatch.setattr(batch, 'build_cloud', lambda args: SimpleNamespace(m=1))
    monkeypatch.setattr(batch, 'build_metric', lambda args, shared: SimpleNamespace(m=1))
    monkeypatch.setattr(batch, 'source_grid', lambda args, cloud: (None, [2., 3.], None))
    def fail(*args, **kwargs):
        raise IndexError('Angular basis too small')
    monkeypatch.setattr(environment_metric_sampling, 'precompute_metric', fail)
    with pytest.raises(IndexError, match='Angular basis too small'):
        batch.main()
    result = json.loads(next(folder.glob('scalar_batch_*.json')).read_text())
    assert result['status'] == 'batch_failed_partial_results_preserved'
    assert result['channels'] == []
    assert result['failure']['stage'] == 'metric_precomputation'


def test_conjugate_only_batch_does_not_solve_direct_channel(tmp_path, monkeypatch):
    folder=tmp_path/'docs'/'environment_reproduction';folder.mkdir(parents=True)
    monkeypatch.setattr(batch,'__file__',str(tmp_path/'src'/'report.py'))
    monkeypatch.setattr(sys,'argv',['report.py','--conjugate-ells','1','--metric-m','2'])
    metric=SimpleNamespace(m=2,_values=SimpleNamespace(cache_info=lambda:SimpleNamespace(hits=0,misses=0)))
    monkeypatch.setattr(batch,'build_cloud',lambda args:SimpleNamespace(m=1))
    monkeypatch.setattr(batch,'build_metric',lambda args,cloud:metric)
    monkeypatch.setattr(batch,'ConjugateMetricMode',lambda base:SimpleNamespace(m=-base.m))
    calls=[]
    def run(args,**kwargs):
        calls.append((args.scalar_ell,args.metric_m))
        return folder/'mode.json',{'flux':{}}
    monkeypatch.setattr(batch,'run',run)
    batch.main()
    assert calls==[(1,-2)]
    result=json.loads(next(folder.glob('scalar_batch_*.json')).read_text())
    assert result['status']=='batch_completed_finite_resolution_not_converged'
    assert result['channels'][0]['scalar_m']==-1


def test_empty_batch_rejected_before_cloud_construction(monkeypatch):
    monkeypatch.setattr(sys,'argv',['report.py'])
    monkeypatch.setattr(batch,'build_cloud',lambda args:pytest.fail('must not build cloud'))
    with pytest.raises(ValueError,match='At least one'):
        batch.main()
