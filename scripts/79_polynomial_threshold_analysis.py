"""Post-hoc reference-guided threshold checks of existing 48 certificates."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import json,hashlib
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent/'run_20261030/polynomial_pinn_certificates';OUT=ROOT.parent/'run_20261030/polynomial_threshold_analysis';OUT.mkdir(parents=True,exist_ok=False)
rows=json.loads((BASE/'records.json').read_text()); tests=[]
for r in rows:
 l,u,rl,ru=map(F,[r['lower'],r['upper'],r['reference_lower'],r['reference_upper']]);mid=(rl+ru)/2
 for distance in (F(1,100),F(1,1000)):
  for direction in (-1,1):
   tau=mid*(1+direction*distance);decision='within' if u<=tau else 'exceeds' if l>tau else 'unknown';reference='within' if ru<=tau else 'exceeds' if rl>tau else 'reference_unknown';ok=decision=='unknown' or reference==decision
   tests.append({'network':r['network'],'source':r['source'],'seed':r['seed'],'degree':r['degree'],'intervals':r['intervals'],'relative_threshold_offset':str(distance),'direction':direction,'threshold':str(tau),'decision':decision,'reference_decision':reference,'compatible':ok})
(OUT/'records.json').write_text(json.dumps(tests,indent=2)+'\n');groups=[]
for d in (2,3):
 for n in (32,64,128,256):
  z=[x for x in tests if x['degree']==d and x['intervals']==n];groups.append({'degree':d,'intervals':n,'total':len(z),'decided':sum(x['decision']!='unknown' for x in z),'unknown':sum(x['decision']=='unknown' for x in z)})
summary={'underlying_certificates':len(rows),'post_hoc_threshold_comparisons':len(tests),'outcomes':dict(Counter(x['decision'] for x in tests)),'compatible':sum(x['compatible'] for x in tests),'failures':[x for x in tests if not x['compatible']],'groups':groups,'source_records_sha256':hashlib.sha256((BASE/'records.json').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'limitation':'thresholds chosen after knowing reference, not blind and no new verification runs'};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k not in ('groups','failures')}),'failures',len(summary['failures']));print('GROUPS',json.dumps(groups))
