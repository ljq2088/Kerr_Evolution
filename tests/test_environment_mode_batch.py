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
