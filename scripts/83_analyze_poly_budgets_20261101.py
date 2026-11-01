"""Reclassify existing decision traces at smaller mesh caps; no new requests."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import json,hashlib
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent/'run_20261101/poly_fixed_tolerance';OUT=ROOT.parent/'run_20261101/poly_budget_analysis';OUT.mkdir(parents=True,exist_ok=False);data=(BASE/'records.json').read_bytes();rows=json.loads(data);budgets=[16,32,64,128];records=[]
for r in rows:
 for cap in budgets:
  trace=[t for t in r['trace'] if t['intervals']<=cap];last=trace[-1] if trace else None
  # Continuing to use an already certified decision is valid at any larger cap.
  records.append({'network':r['network'],'source':r['source'],'degree':r['degree'],'tolerance':r['tolerance'],'max_intervals':cap,'decision':last['decision'] if last else 'no_result','used_intervals':last['intervals'] if last else None})
groups=[]
for cap in budgets:
 for degree in (2,3):
  group=[x for x in records if x['max_intervals']==cap and x['degree']==degree];groups.append({'max_intervals':cap,'degree':degree,'total':len(group),'outcomes':dict(Counter(x['decision'] for x in group))})
summary={'underlying_requests':len(rows),'derived_budget_classifications':len(records),'groups':groups,'data_sha256':hashlib.sha256(data).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'limitation':'derived retrospectively from existing traces; not new independent runs'};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');(OUT/'records.json').write_text(json.dumps(records,indent=2)+'\n');print(json.dumps(summary,indent=2))
