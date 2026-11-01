"""Validate completed 106/107 research outputs in a separate process.

Checks storage consistency and request-key coverage, not mathematical validity.
Only the exact record schemas emitted by scripts 106 and 107 are supported.
"""
import argparse,hashlib,json,sys
from pathlib import Path

def verify(directory, fixtures=None):
 errors=[]
 try:
  raw=(directory/'records.json').read_bytes()
  rows=json.loads(raw)
  summary=json.loads((directory/'summary.json').read_text())
  if not isinstance(rows,list) or not isinstance(summary,dict):raise ValueError('records must be an array and summary an object')
  digest=hashlib.sha256(raw).hexdigest()
  if digest!=summary['records_sha256']:errors.append('records SHA-256 mismatch')
  if len(rows)!=summary['records'] or len(rows)!=summary['expected']:errors.append('record count mismatch')
  if 'fixtures' in summary and 'points' in summary and 'orders' in summary:
   kind='derivative'; keys=[(x['fixture_sha256'],x['network'],x['x'],x['order']) for x in rows]
   declared={f['sha256']:f['networks'] for f in summary['fixtures']}
   if len(declared)!=len(summary['fixtures']):errors.append('duplicate fixture hash')
   for sha,num in declared.items():
    group=[x for x in rows if x['fixture_sha256']==sha]
    labels={x['network'] for x in group}
    expected={(sha,n,p,o) for n in labels for p in summary['points'] for o in summary['orders']}
    if len(labels)!=num or set(k for k in keys if k[0]==sha)!=expected:errors.append('fixture request coverage mismatch')
   if set(k[0] for k in keys)!=set(declared):errors.append('unexpected fixture')
   if summary['stable']!=sum(x['mp_60_vs_90_stable'] is True for x in rows):errors.append('stability total mismatch')
   if summary['contained']!=sum(x['mp90_contained'] is True for x in rows):errors.append('containment total mismatch')
  elif (directory/'config.json').exists():
   kind='certificate';config=json.loads((directory/'config.json').read_text());keys=[(x['network'],x['degree']) for x in rows];labels={x['network'] for x in rows};expected={(n,d) for n in labels for d in config['degrees']}
   if len(labels)!=config['networks'] or set(keys)!=expected:errors.append('network/degree coverage mismatch')
   if any(x['intervals']!=config['intervals'] for x in rows):errors.append('unexpected mesh')
   if summary['contained']!=sum(x['reference_contained'] is True for x in rows):errors.append('containment total mismatch')
   if summary['overlap']!=sum(x['reference_overlap'] is True for x in rows):errors.append('overlap total mismatch')
  else:raise ValueError('unrecognized result schema')
  if len(keys)!=len(set(keys)):errors.append('duplicate request key')
  if fixtures:
   actual={}
   for path in fixtures:
    raw_fixture=path.read_bytes();sha=hashlib.sha256(raw_fixture).hexdigest();nets=json.loads(raw_fixture)['networks'];labels=[n['label'] for n in nets]
    if len(set(labels))!=len(labels):errors.append('duplicate fixture network label')
    if sha in actual:errors.append('duplicate supplied fixture')
    actual[sha]=set(labels)
   if kind=='certificate':
    if set(actual)!={config['fixture_sha256']}:errors.append('fixture SHA-256 mismatch')
    else:
     manifest={(n,d) for n in actual[config['fixture_sha256']] for d in config['degrees']}
     if set(keys)!=manifest:errors.append('network keys differ from supplied fixture')
   else:
    if set(actual)!=set(declared):errors.append('fixture SHA-256 mismatch')
    else:
     manifest={(sha,n,p,o) for sha,labels in actual.items() for n in labels for p in summary['points'] for o in summary['orders']}
     if set(keys)!=manifest:errors.append('network keys differ from supplied fixture')
  if summary.get('failures'):errors.append('experiment reports failures')
  return {'ok':not errors,'directory':str(directory),'schema':kind,'records':len(rows),'sha256':digest,'errors':errors,'scope':'record consistency and coverage, not mathematical validity'}
 except (OSError,ValueError,KeyError,TypeError,json.JSONDecodeError) as e:
  return {'ok':False,'directory':str(directory),'errors':[f'{type(e).__name__}: {e}'],'scope':'record consistency and coverage, not mathematical validity'}

def main():
 p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('--fixtures',type=Path,nargs='+');a=p.parse_args();r=verify(a.directory,a.fixtures);print(json.dumps(r,ensure_ascii=False));return 0 if r['ok'] else 1
if __name__=='__main__':sys.exit(main())
