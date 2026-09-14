#!/usr/bin/env python3
import csv, json, re, sys
from collections import Counter
from pathlib import Path

root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parents[1]
expected={'PRJEB32839':(150,75),'PRJEB14695':(138,23),'PRJNA597909':(20,20),'PRJNA602528':(10,10)}
errors=[]
for project,(nr,ns) in expected.items():
    study=root/'studies'/project
    runs=list(csv.DictReader((study/'runs.tsv').open(encoding='utf-8'),delimiter='	'))
    samples=list(csv.DictReader((study/'samples.tsv').open(encoding='utf-8'),delimiter='	'))
    metadata=list(csv.DictReader((study/'metadata.csv').open(encoding='utf-8')))
    if (len(runs),len(samples),len(metadata)) != (nr,ns,nr): errors.append(f'{project}: counts')
    if len({r['run_accession'] for r in runs}) != nr: errors.append(f'{project}: duplicate runs')
    if {r['run_accession'] for r in runs} != {r['run_accession'] for r in metadata}: errors.append(f'{project}: run mismatch')
    if {s['sample_id'] for s in samples} != {r['sample_id'] for r in runs}: errors.append(f'{project}: sample mismatch')
    counts=Counter(r['sample_id'] for r in runs)
    for s in samples:
        if counts[s['sample_id']] != int(s['technical_runs']): errors.append(f"{project}: technical count {s['sample_id']}")
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]*',s['sample_id']): errors.append(f"{project}: unsafe id {s['sample_id']}")
    if any(r['library_layout']!='PAIRED' for r in runs): errors.append(f'{project}: non-paired run')
    if project != 'PRJNA602528':
        spec=json.loads((study/'de_spec.json').read_text(encoding='utf-8'))
        levels=Counter(s['condition'] for s in samples)
        for contrast in spec['contrasts']:
            for level in (contrast['numerator'],contrast['denominator']):
                if levels[level] < spec['parameters']['min_replicates']: errors.append(f'{project}: under-replicated {level}')
    else:
        blocked=json.loads((study/'de_spec.blocked.json').read_text(encoding='utf-8'))
        if blocked.get('status')!='BLOCKED': errors.append(f'{project}: DE block missing')
print(json.dumps({'status':'PASS' if not errors else 'FAIL','errors':errors},indent=2))
raise SystemExit(1 if errors else 0)
