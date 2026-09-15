"""Fetch the pinned, small official Mathematica paclets into ignored outputs.
No global package installation. Existing files must match; nothing is overwritten.
"""
from pathlib import Path
import hashlib
import io
import json
import tarfile
import urllib.request

ROOT = Path(__file__).resolve().parents[3]
MANIFEST = Path(__file__).with_name('official_bhpt_manifest.json')
DEST = ROOT / 'outputs/paper_metric_reference/wolfram_dependencies'

def main():
    DEST.mkdir(parents=True, exist_ok=True)
    for repository in json.loads(MANIFEST.read_text())['repositories']:
        name = repository['repository'].split('/')[-1]
        archive = DEST / (name + '_' + repository['commit'] + '.tar.gz')
        if archive.exists():
            data = archive.read_bytes()
        else:
            request = urllib.request.Request(repository['archive_url'], headers={'User-Agent':'BHPT-independent-audit'})
            with urllib.request.urlopen(request, timeout=30) as response:
                data = response.read(4_000_001)
        if len(data) != repository['archive_bytes'] or hashlib.sha256(data).hexdigest() != repository['archive_sha256']:
            raise ValueError('Pinned archive checksum mismatch: ' + name)
        if not archive.exists():
            archive.write_bytes(data)
        target = DEST / name
        target.mkdir(exist_ok=True)
        with tarfile.open(fileobj=io.BytesIO(data)) as tar:
            for member in tar:
                relative = Path(*Path(member.name).parts[1:])
                if not relative.parts:
                    continue
                path = (target / relative).resolve()
                if target.resolve() not in path.parents:
                    raise ValueError('Unsafe archive member')
                if member.isdir():
                    path.mkdir(parents=True, exist_ok=True)
                elif member.isfile():
                    payload = tar.extractfile(member).read()
                    path.parent.mkdir(parents=True, exist_ok=True)
                    if path.exists():
                        if path.read_bytes() != payload:
                            raise ValueError('Existing file differs: ' + str(path))
                    else:
                        path.write_bytes(payload)
                else:
                    raise ValueError('Unexpected archive member type')
        print(name, repository['commit'], 'verified')

if __name__ == '__main__':
    main()
