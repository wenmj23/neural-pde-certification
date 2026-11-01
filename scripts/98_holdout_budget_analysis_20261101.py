"""Retrospective mesh-cap outcomes from prelocked holdout traces; no new runs."""
from pathlib import Path
from collections import Counter
import json,hashlib
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent/'run_20261101/holdout_benchmark';OUT=ROOT.parent/'run_20261101/holdout_budget_analysis';OUT.mkdir(parents=True,exist_ok=False);data=(BASE/'records.json').read_bytes();rows=json.loads(data);checks=json.loads((BASE/'checks.json').read_text());ref={(x['network'],x['tolerance'],x['method'],x['repeat']):x['reference_class'] for x in checks};out=[]
for r in rows:
 for cap in (16,32,64,128):
  eligible=[x for x in r['trace'] if x['intervals']<=cap];last=eligible[-1] if eligible else None;reference=ref[(r['network'],r['tolerance'],r['method'],r['repeat'])];status='within' if last and last['within'] else 'unknown';out.append({'network':r['network'],'tolerance':r['tolerance'],'method':r['method'],'repeat':r['repeat'],'cap':cap,'used_intervals':last['intervals'] if last else None,'decision':status,'reference_class':reference})
groups=[]
for tau in ('1/1000','1/100000'):
 for cap in (16,32,64,128):
  for method in ('global','spatial','signed_quadratic','signed_cubic'):
   z=[r for r in out if r['tolerance']==tau and r['cap']==cap and r['method']==method and r['reference_class']=='within'];groups.append({'tolerance':tau,'cap':cap,'method':method,'eligible_requests':len(z),'within_requests':sum(r['decision']=='within' for r in z),'unknown_requests':sum(r['decision']=='unknown' for r in z)})
summary={'source_requests':len(rows),'derived_classifications':len(out),'groups':groups,'data_sha256':hashlib.sha256(data).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'note':'reclassified existing traces; repeated counts are not independent networks, and an early successful certificate remains valid at larger caps'};(OUT/'records.json').write_text(json.dumps(out,indent=2)+'\n');(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({'source_requests':len(rows),'derived_classifications':len(out),'groups':groups},indent=2))
