"""Increase budgets only for the unresolved 0.01% cases from script 56; post-hoc follow-up."""
import json, subprocess, sys, time, hashlib
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; PRE=ROOT.parent/'run_20261030'/'budget_curves'; OUT=ROOT.parent/'run_20261030'/'budget_followup'; OUT.mkdir(parents=True,exist_ok=False)
prior=json.loads((PRE/'records.json').read_text()); todo=[r for r in prior if r['relative_distance']=='1/10000' and r['cap']==256 and r['decision']=='undetermined']; rows=[]
def f(x):return F(int(x['numerator']),int(x['denominator']))
for r in todo:
 for cap in (512,1024):
  tag=f"{r['case']}_side{r['direction']:+d}_degree{r['degree']}_cap{cap}"; dst=OUT/tag
  cmd=[sys.executable,str(ROOT/'scripts/certify.py'),'--case',str(ROOT/f"configs/examples/{r['case']}.json"),'--out',str(dst),'--degree',str(r['degree']),'--mesh','uniform','--goal','decision','--error-tolerance',r['threshold'],'--max-intervals',str(cap),'--max-attempts','9','--time-limit-seconds','30']
  t=time.perf_counter(); proc=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,timeout=60); path=dst/'result.json'; data=json.loads(path.read_text()) if path.exists() else None; enc=data['enclosure'] if data else None; dec=data['accuracy_decision'] if data else 'missing'
  row={'case':r['case'],'direction':r['direction'],'degree':r['degree'],'cap':cap,'threshold':r['threshold'],'decision':dec,'stop_reason':data['stop_reason'] if data else 'missing','intervals':enc['intervals'] if enc else None,'compatible':dec in ('undetermined',r['expected_from_reference']),'reference_contained':bool(enc and f(enc['lower'])<=F(r['reference_lower']) and F(r['reference_upper'])<=f(enc['upper'])),'returncode':proc.returncode,'wall_seconds_external':time.perf_counter()-t,'stdout':proc.stdout,'stderr':proc.stderr}
  rows.append(row); (OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n');print(tag,dec,row['intervals'],flush=True)
summary={'inputs':len(todo),'requests':len(rows),'compatible':sum(r['compatible'] for r in rows),'reference_contained':sum(r['reference_contained'] for r in rows),'determined':sum(r['decision']!='undetermined' for r in rows),'errors':[r for r in rows if not r['compatible'] or not r['reference_contained'] or r['returncode'] not in (0,2)],'prior_records_sha256':hashlib.sha256((PRE/'records.json').read_bytes()).hexdigest()}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',summary,flush=True)
