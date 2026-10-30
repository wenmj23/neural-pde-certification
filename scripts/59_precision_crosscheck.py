"""Cross-check three allowed precisions; same core, not independent audit."""
import json,hashlib,subprocess,sys,time
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
ROOT=Path(__file__).resolve().parents[1]; BASE=ROOT.parent/'run_20261030/frozen_portfolio'; OUT=ROOT.parent/'run_20261030/precision_crosscheck'; OUT.mkdir(parents=True,exist_ok=False)
names=['depth1_k1_seed7','depth1_k5_seed37','depth2_k1_seed17','depth3_k5_seed7','depth1_k2_seed49','depth3_k5_seed49']; refs=json.loads((ROOT.parent/'run_20261030/public_crosscheck/published_references.json').read_text());rows=[];checks=[]
def fr(d):return F(int(d['numerator']),int(d['denominator']))
for name in names:
 case=json.loads((BASE/(name+'.case.json')).read_text()); pair=[]
 for bits in (80,103,160):
  case['verification']['degree']=3;case['verification']['precision_bits']=bits
  file=OUT/f'{name}__p{bits}.case.json';file.write_text(json.dumps(case,indent=2)+'\n'); dest=OUT/f'{name}__p{bits}';cmd=[sys.executable,str(ROOT/'scripts/certify.py'),'--case',str(file),'--out',str(dest)]
  t=time.perf_counter();p=subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True,timeout=60);rp=dest/'result.json';r=json.loads(rp.read_text()) if rp.exists() else None;e=r['enclosure'] if r else None;ref=refs.get(name)
  item={'name':name,'bits':bits,'returncode':p.returncode,'goal_met':r['goal_met'] if r else False,'stop_reason':r['stop_reason'] if r else 'no_result','intervals':e['intervals'] if e else None,'lower':str(fr(e['lower'])) if e else None,'upper':str(fr(e['upper'])) if e else None,'reference_contained':bool(fr(e['lower'])<=fr(ref['lower']) and fr(ref['upper'])<=fr(e['upper'])) if e and ref else None,'wall_seconds_external':time.perf_counter()-t,'stdout':p.stdout,'stderr':p.stderr};rows.append(item);pair.append(item);(OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n');print(name,bits,item['stop_reason'],item['intervals'],flush=True)
 for a,b in combinations(pair,2):
  overlap=(max(F(a['lower']),F(b['lower']))<=min(F(a['upper']),F(b['upper']))) if a['lower'] and b['lower'] else None
  checks.append({'name':name,'precision_pair':[a['bits'],b['bits']],'overlap':overlap})
summary={'requests':len(rows),'networks':len(names),'goal_met':sum(r['goal_met'] for r in rows),'reference_checked':sum(r['reference_contained'] is not None for r in rows),'reference_contained':sum(r['reference_contained'] is True for r in rows),'overlap_checked':len(checks),'overlap_passed':sum(c['overlap'] is True for c in checks),'unexpected':[r for r in rows if r['returncode'] not in (0,2)],'checks':checks,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k not in ('checks','unexpected')}),'unexpected',len(summary['unexpected']),flush=True)
