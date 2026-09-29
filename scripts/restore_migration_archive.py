"""Verify a migration archive; optionally restore missing, byte-identical files.

Different existing files are never overwritten. Extraction uses only regular
files and rejects traversal, links and symlinked parents. The first pass verifies
all data before any destination writes; a failure never silently replaces work.
"""
import argparse,hashlib,json,tarfile
from pathlib import Path,PurePosixPath

def sha(stream):
    h=hashlib.sha256()
    for chunk in iter(lambda:stream.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def target(root,name):
    rel=PurePosixPath(name)
    if rel.is_absolute() or '..' in rel.parts or '\\' in name or ':' in name: raise ValueError('Unsafe member: '+name)
    path=root.joinpath(*rel.parts)
    if not path.resolve().is_relative_to(root): raise ValueError('Escaping member: '+name)
    for part in [path,*path.parents]:
        if part==root:break
        if part.is_symlink():raise ValueError('Symlink destination: '+str(part))
    return path

def main():
    ap=argparse.ArgumentParser();ap.add_argument('archive',type=Path);ap.add_argument('--destination',type=Path)
    a=ap.parse_args();root=a.destination.resolve() if a.destination else None
    with tarfile.open(a.archive,'r:gz') as tar:
        manifest=json.load(tar.extractfile('MIGRATION_FILES.json'));expected={r['path']:r for r in manifest['files']};seen=set()
        for member in tar:
            if member.name=='MIGRATION_FILES.json':continue
            if member.name in seen or member.name not in expected or not member.isfile():raise ValueError('Unexpected/duplicate member: '+member.name)
            r=expected[member.name];seen.add(member.name)
            if member.size!=r['size'] or sha(tar.extractfile(member))!=r['sha256']:raise ValueError('Checksum failed: '+member.name)
            if root:
                out=target(root,member.name)
                if out.exists():
                    if not out.is_file():raise ValueError('Existing non-file: '+str(out))
                    with out.open('rb') as f:
                        if sha(f)!=r['sha256']:raise ValueError('Different existing file; use a fresh destination: '+str(out))
        if seen!=set(expected):raise ValueError('Missing members')
    if root:
        root.mkdir(parents=True,exist_ok=True)
        with tarfile.open(a.archive,'r:gz') as tar:
            for member in tar:
                if member.name=='MIGRATION_FILES.json':continue
                out=target(root,member.name)
                if out.exists():continue
                out.parent.mkdir(parents=True,exist_ok=True)
                with out.open('xb') as f,tar.extractfile(member) as source:
                    for chunk in iter(lambda:source.read(1024*1024),b''):f.write(chunk)
                out.chmod(member.mode & 0o777)
    print(json.dumps({'verified_files':len(seen),'restored_to':str(root) if root else None}))
if __name__=='__main__':main()
