"""Restore archived evidence to the generated paths expected by old experiments.

Copies verified bytes without rerunning training or experiments. Existing files
with different contents are rejected. Use --destination-root for an isolated
restore; the default is the project root. No network or credentials are used.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

ROOT=Path(__file__).resolve().parents[1]
MANIFEST=ROOT/'results/2026-10-28/research_runs/archive_manifest.json'


def restore(destination_root):
    destination_root=Path(destination_root).resolve()
    manifest=json.loads(MANIFEST.read_text())
    plan=[]
    for entry in manifest['entries']:
        source=(ROOT/entry['archive_path']).resolve()
        target=(destination_root/entry['generated_path']).resolve()
        if not source.is_relative_to(ROOT/'results/2026-10-28/research_runs'):
            raise ValueError('Archive path escapes evidence directory')
        if not target.is_relative_to(destination_root/'results/generated'):
            raise ValueError('Restore path escapes generated directory')
        digest=hashlib.sha256(source.read_bytes()).hexdigest()
        if digest!=entry['sha256']:raise ValueError('Archive checksum mismatch: '+entry['archive_path'])
        if target.exists() and hashlib.sha256(target.read_bytes()).hexdigest()!=digest:
            raise ValueError('Existing generated file differs; preserve it separately before restoring: '+str(target))
        plan.append((source,target))
    copied=0
    for source,target in plan:
        if not target.exists():
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(source,target);copied+=1
    return dict(verified=len(plan),copied=copied,already_present=len(plan)-copied)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--destination-root',type=Path,default=ROOT)
    args=parser.parse_args()
    print(json.dumps(restore(args.destination_root),indent=2))
