"""Post-hoc threshold stress test; not blind evaluation."""
import json, os, subprocess, sys
from pathlib import Path
from fractions import Fraction
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT.parent/'run_20261030'/'public_crosscheck'
cases=['shallow_sine','deep_sine','holdout_shallow_k2','holdout_deep_k5']; rows=[]
def f(x):return Fraction(int(x['numerator']),int(x['denominator']))
for c in cases:
 b=json.loads((OUT/(c+'__budget_16')/'result.json').read_text())['enclosure']; lo,hi=f(b['lower']),f(b['upper']); tau=(lo+hi)/2
 dest=OUT/(c+'__decision_straddle_16'); cmd=[sys.executable,str(ROOT/'scripts/certify.py'),'--case',str(ROOT/f'configs/examples/{c}.json'),'--out',str(dest),'--degree','3','--goal','decision','--max-intervals','16','--max-attempts','1','--error-tolerance',str(tau)]
 p=subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True,timeout=90); d=json.loads((dest/'result.json').read_text()) if (dest/'result.json').exists() else None
 ok=bool(d and not d['goal_met'] and d['accuracy_decision']=='undetermined' and p.returncode==2)
 rows.append({'case':c,'threshold':str(tau),'selection':'midpoint of previous 16-cell enclosure; post-hoc','returncode':p.returncode,'passed':ok,'stdout':p.stdout,'stderr':p.stderr})
 print(c,'unknown',ok,flush=True)
(OUT/'straddle_checks.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2)+'\n')
