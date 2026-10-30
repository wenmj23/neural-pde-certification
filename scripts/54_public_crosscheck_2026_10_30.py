"""Finite cross-check of public CLI; not an independent software audit."""
import json, os, platform, subprocess, sys, time
from fractions import Fraction
from pathlib import Path
from itertools import combinations
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT.parent/'run_20261030'/'public_crosscheck'
OUT.mkdir(parents=True, exist_ok=False)
cases=['shallow_sine','deep_sine','holdout_shallow_k2','holdout_deep_k5']
records=[]
def run(case, tag, extra):
    dest=OUT/f'{case}__{tag}'
    cmd=[sys.executable,str(ROOT/'scripts/certify.py'),'--case',str(ROOT/f'configs/examples/{case}.json'),'--out',str(dest),*extra]
    tick=time.perf_counter(); p=subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True,timeout=90)
    data=json.loads((dest/'result.json').read_text()) if (dest/'result.json').exists() else None
    row={'case':case,'tag':tag,'arguments':extra,'returncode':p.returncode,'elapsed_total_seconds':time.perf_counter()-tick,'stdout':p.stdout,'stderr':p.stderr,'result_path':str(dest/'result.json') if data else None}
    records.append(row); (OUT/'runs.json').write_text(json.dumps(records,indent=2,ensure_ascii=False)+'\n')
    print(case,tag,'exit',p.returncode,'stop',data['stop_reason'] if data else 'no_result',flush=True)
    return data

def frac(x): return Fraction(int(x['numerator']),int(x['denominator']))
checks=[]
for case in cases:
    results=[]
    for degree in (2,3):
        for mesh in ('uniform','adaptive'):
            d=run(case,f'd{degree}_{mesh}',['--degree',str(degree),'--mesh',mesh])
            if d and d['enclosure']: results.append((degree,mesh,d))
            checks.append({'case':case,'check':f'd{degree}_{mesh}_goal_met','passed':bool(d and d['goal_met'])})
    for (d1,m1,a),(d2,m2,b) in combinations(results,2):
        a,b=a['enclosure'],b['enclosure']
        passed=max(frac(a['lower']),frac(b['lower']))<=min(frac(a['upper']),frac(b['upper']))
        checks.append({'case':case,'check':f'overlap_d{d1}_{m1}_d{d2}_{m2}','passed':passed})
    d=run(case,'budget_16',['--degree','3','--goal','enclosure','--max-intervals','16','--max-attempts','1','--gap-ratio','1.000000000001'])
    checks.append({'case':case,'check':'limited_budget_not_success','passed':bool(d and not d['goal_met'] and d['stop_reason'] in ('mesh_limit','attempt_limit'))})
    for tolerance in ('100','1e-12'):
        d=run(case,'decision_'+('large' if tolerance=='100' else 'small'),['--degree','3','--goal','decision','--error-tolerance',tolerance])
        expected='within_tolerance' if tolerance=='100' else 'exceeds_tolerance'
        checks.append({'case':case,'check':'decision_'+tolerance,'passed':bool(d and d['goal_met'] and d['accuracy_decision']==expected)})
summary={'repository_commit':'1948ea20ae6d88395dc417c97581a78a427484da','python':platform.python_version(),'case_count':len(cases),'request_count':len(records),'checks':checks,'check_count':len(checks),'passed_count':sum(c['passed'] for c in checks),'failed_checks':[c for c in checks if not c['passed']]}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='checks'},ensure_ascii=False),flush=True)
