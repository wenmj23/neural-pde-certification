"""Check persisted request count, uniqueness, decisions and reference comparisons."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import argparse,json,sys
p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('--output',type=Path);args=p.parse_args();root=args.directory
rows=json.loads((root/'records.json').read_text());summary=json.loads((root/'summary.json').read_text());config=json.loads((root/'config.json').read_text());keys=[(r['network'],r['degree'],r['tolerance']) for r in rows];errors=[]
if len(rows)!=72:errors.append(f'expected 72 records, got {len(rows)}')
if len(set(keys))!=len(keys):errors.append('duplicate request keys')
if summary['requests']!=len(rows):errors.append(f"summary says {summary['requests']} but records contain {len(rows)}")
if dict(Counter(r['decision'] for r in rows))!=summary['outcomes']:errors.append('outcome counts disagree')
for r in rows:
 tau=F(r['tolerance']);trace=r['trace']
 if not trace:errors.append('empty trace: '+str((r['network'],r['degree'],r['tolerance'])));continue
 for t in trace:
  lo,hi=F(t['lower']),F(t['upper']);expected='within' if hi<=tau else 'exceeds' if lo>tau else 'unknown'
  if lo<0 or hi<lo or t['decision']!=expected:errors.append('invalid trace: '+str((r['network'],r['degree'],r['tolerance'])))
 if r['decision']!=trace[-1]['decision']:errors.append('final differs from trace: '+str((r['network'],r['degree'],r['tolerance'])))
 if not r['compatible']:errors.append('reference disagreement: '+str((r['network'],r['degree'],r['tolerance'])))
result={'ok':not errors,'record_count':len(rows),'summary_count':summary['requests'],'unique_count':len(set(keys)),'errors':errors}
if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,ensure_ascii=False));sys.exit(0 if not errors else 1)
