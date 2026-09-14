#!/usr/bin/env python3
import argparse, csv, hashlib, json
from pathlib import Path

def md5(path):
    h=hashlib.md5()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()

p=argparse.ArgumentParser(); p.add_argument('project'); p.add_argument('--package-root',type=Path,required=True); p.add_argument('--data-root',type=Path,required=True); a=p.parse_args()
manifest=a.package_root/'studies'/a.project/'runs.tsv'
rows=list(csv.DictReader(manifest.open(encoding='utf-8'),delimiter='	'))
errors=[]
for row in rows:
    for mate in ('1','2'):
        path=a.data_root/a.project/'fastq_ftp'/f"{row['run_accession']}_{mate}.fastq.gz"
        if not path.is_file(): errors.append(f'missing: {path}'); continue
        expected_bytes=int(row[f'fastq_{mate}_bytes'])
        if path.stat().st_size != expected_bytes: errors.append(f'bytes: {path}')
        if md5(path) != row[f'fastq_{mate}_md5']: errors.append(f'md5: {path}')
print(json.dumps({'project':a.project,'runs':len(rows),'status':'PASS' if not errors else 'FAIL','errors':errors},indent=2))
raise SystemExit(1 if errors else 0)
