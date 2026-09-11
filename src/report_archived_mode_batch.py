"""Run an explicitly pinned historical solver against its own metric cache.

Cached tensor files and metadata are never rewritten. All numerical modules
come from the requested Git commit; only report output paths target this repo.
"""
import argparse
import hashlib
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-commit', required=True)
    parser.add_argument('--expected-source-hash', required=True)
    args, remaining = parser.parse_known_args()
    if remaining[:1] == ['--']:
        remaining = remaining[1:]
    repo = Path(__file__).resolve().parents[1]
    revision = subprocess.check_output(['git', 'rev-parse', '--verify',
        args.source_commit+'^{commit}'], cwd=repo, text=True).strip()
    source = repo / 'outputs' / 'static_recovery_source' / revision / 'src'
    names = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only',
        revision, 'src'], cwd=repo, text=True).splitlines()
    for name in names:
        if not name.endswith('.py'):
            continue
        relative = Path(name).relative_to('src')
        if '..' in relative.parts:
            raise ValueError('Unsafe archived source path')
        target = source / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(subprocess.check_output(['git', 'show', revision+':'+name], cwd=repo))
    paths = sorted(set(source.glob('lorenz_*.py')) | {source / name for name in (
        'environment_lorenz_mode.py', 'environment_static_lorenz.py',
        'environment_radial.py', 'environment_source.py', 'environment_cloud.py')})
    digest = hashlib.sha256()
    for path in paths:
        digest.update(path.name.encode()); digest.update(path.read_bytes())
    if digest.hexdigest() != args.expected_source_hash:
        raise ValueError('Historical numerical source hash does not match the requested cache version')
    sys.path.insert(0, str(source))
    import report_environment_mode_batch as batch
    import report_environment_forced_mode as forced
    # Numerical modules keep their historical __file__, including cache hashing.
    batch.__file__ = str(repo / 'src' / 'report_environment_mode_batch.py')
    forced.__file__ = str(repo / 'src' / 'report_environment_forced_mode.py')
    sys.argv = ['report_environment_mode_batch.py', *remaining]
    print(f'Historical solver commit: {revision}; numerical source hash: {digest.hexdigest()}', flush=True)
    batch.main()


if __name__ == '__main__':
    main()
