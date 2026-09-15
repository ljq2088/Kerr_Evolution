#!/usr/bin/env python3
"""Rebuild Sam audit inputs without redistributing the author's source.

Requires the fixed public Sam .m and the adjacent versioned audit manifests.
This generates files only; it never starts Mathematica. For i7, provide the
completed i6 matching_state.mx. A fresh checkpoint can have a different binary
hash, so that hash is inserted into our driver and recorded in the build record.
"""
import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
MANIFESTS = HERE.parent
SOURCE_SHA256 = '124fc609fb36f46d4f28378ac040655c8ae7d3cff9346ec62134c358ce223ca9'
SOURCE_COMMIT = '5c1b42793ff893fd0c65a5b48ceecc02d34c5e35'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def manifest(suffix):
    return json.loads((MANIFESTS/f'sam_actual_{suffix}_driver_manifest_20260915.json').read_text())


def build(args):
    original = args.source.read_bytes()
    if sha(original) != SOURCE_SHA256:
        raise ValueError('Source is not the fixed Sam commit/hash; obtain the documented public archive')
    i6 = manifest('i6_h5')
    selected = manifest(args.variant)
    lines = original.decode('utf-8').splitlines()
    files = {}
    for name, meta in i6['phase_sources'].items():
        parts = [f'(* Original source lines {lo}-{hi}, unchanged. *)\n'+'\n'.join(lines[lo-1:hi])
                 for lo,hi in meta['original_line_spans']]
        text = '\n\n'.join(parts)+'\n'
        if name == '01_matching.m':
            old = 'a0=SetPrecision[Round[GetParam["a"]*iround]/iround,prec];'
            new = 'a0=SetPrecision[GetParam["a"],prec]; (* Audit driver: preserve supplied exact a, no 3-digit rounding. *)'
            if text.count(old) != 1:
                raise ValueError('Unexpected original a0 initialization')
            text = text.replace(old,new)
        elif name == '02_radial.m':
            before,after = {},{}
            for change in i6['changes_in_this_run']['logging_insertions']:
                if 'before_source_line' in change:
                    before.setdefault(change['before_source_line'],[]).append(change['added'])
                else:
                    after.setdefault(change['after_source_line'],[]).append(change['added'])
            modified=[]
            # First generated line is our label; next line is original line954.
            for line_no,line in enumerate(text.splitlines(),start=953):
                modified.extend(before.get(line_no,[]));modified.append(line);modified.extend(after.get(line_no,[]))
            text='\n'.join(modified)+'\n'
        data=text.encode('utf-8')
        if sha(data) != meta['sha256']:
            raise ValueError('Regenerated source block hash differs: '+name)
        files[name]=data
    config=''.join(key+'\t'+value+'\n' for key,value in selected['config'].items())
    files['config/configaudit.txt']=config.encode('utf-8')
    template=(HERE/f'run_sam_{args.variant}.wls').read_bytes()
    if sha(template) != selected['driver_sha256']:
        raise ValueError('Versioned orchestration template differs from run manifest')
    if args.variant == 'i7_h6':
        if args.checkpoint is None:
            raise ValueError('i7_h6 requires --checkpoint from a completed i6_h5 run')
        checkpoint=args.checkpoint.read_bytes()
        old=selected['restored_matching_state']['sha256']
        new=sha(checkpoint)
        marker=f'AuthorAudit`expectedCheckpointSHA="{old}";'.encode()
        if template.count(marker)!=1:
            raise ValueError('Unexpected checkpoint verification marker')
        template=template.replace(marker,f'AuthorAudit`expectedCheckpointSHA="{new}";'.encode())
        files['matching_state.mx']=checkpoint
    files['run_sam.wls']=template
    args.output.mkdir(parents=True,exist_ok=True)
    for name,data in files.items():
        target=args.output/name
        target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes(data)
    (args.output/'data').mkdir(exist_ok=True)
    record={'status':'generated_only_not_executed','variant':args.variant,'source_commit':SOURCE_COMMIT,
            'source_sha256':SOURCE_SHA256,'source_path':str(args.source.resolve()),
            'source_line_spans':{n:m['original_line_spans'] for n,m in i6['phase_sources'].items()},
            'config':selected['config'],'files_sha256':{n:sha(data) for n,data in files.items()},
            'driver_template_sha256':sha((HERE/f'run_sam_{args.variant}.wls').read_bytes()),
            'generator_sha256':sha(Path(__file__).read_bytes()),
            'dependency_setup':'Place the documented bhpt_paths.wl one directory above the output folder; run_sam.wls loads it.'}
    (args.output/'build_manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'output':str(args.output),'variant':args.variant,'source_blocks_verified':True}))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=ROOT/'outputs/paper_metric_reference/srd24_KerrLorenzCirc/metric_reconstruction_calc_radiative.m')
    p.add_argument('--variant',choices=['i6_h5','i7_h6'],required=True)
    p.add_argument('--checkpoint',type=Path)
    p.add_argument('--output',type=Path,required=True)
    build(p.parse_args())


if __name__=='__main__':main()
