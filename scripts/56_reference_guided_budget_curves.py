"""Reference-guided decision budget experiment; not a blind holdout or proof."""
import json, subprocess, sys, time, hashlib, platform
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT.parent/'run_20261030'/'public_crosscheck'
OUT=ROOT.parent/'run_20261030'/'budget_curves'
OUT.mkdir(parents=True,exist_ok=False)
refs=json.loads((BASE/'published_references.json').read_text())
cases={'holdout_shallow_k2':'depth1_k2_seed49','holdout_deep_k5':'depth3_k5_seed49'}
def frac(d): return F(int(d['numerator']),int(d['denominator']))
def save(data): (OUT/'records.json').write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
rows=[]
for case,label in cases.items():
 ref=refs[label]; rl,ru=frac(ref['lower']),frac(ref['upper']); mid=(rl+ru)/2
 for distance in (F(1,100),F(1,1000),F(1,10000)):
  for direction in (-1,1):
   tau=mid*(1+direction*distance)
   for degree in (2,3):
    for cap in (64,256):
     tag=f'{case}_dist{distance.denominator}_side{direction:+d}_degree{degree}_cap{cap}'
     dst=OUT/tag; cmd=[sys.executable,str(ROOT/'scripts/certify.py'),'--case',str(ROOT/f'configs/examples/{case}.json'),'--out',str(dst),'--degree',str(degree),'--mesh','uniform','--goal','decision','--error-tolerance',str(tau),'--max-intervals',str(cap),'--max-attempts','9','--time-limit-seconds','30']
     t=time.perf_counter(); proc=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,timeout=60)
     path=dst/'result.json'; res=json.loads(path.read_text()) if path.exists() else None
     enc=res['enclosure'] if res else None
     dec=res['accuracy_decision'] if res else 'missing'
     correct=('exceeds_tolerance' if ru>rl and tau<rl else 'within_tolerance' if tau>=ru else 'reference_overlaps_threshold')
     compatible=(dec=='undetermined' or dec==correct) if res else False
     contained=(frac(enc['lower'])<=rl and ru<=frac(enc['upper'])) if enc else False
     row={'case':case,'label':label,'relative_distance':str(distance),'direction':direction,'degree':degree,'cap':cap,'threshold':str(tau),'reference_lower':str(rl),'reference_upper':str(ru),'expected_from_reference':correct,'decision':dec,'decision_compatible':compatible,'reference_contained':contained,'goal_met':res['goal_met'] if res else False,'stop_reason':res['stop_reason'] if res else 'missing','intervals':enc['intervals'] if enc else None,'returncode':proc.returncode,'wall_seconds_external':time.perf_counter()-t,'stdout':proc.stdout,'stderr':proc.stderr,'output_directory':str(dst)}
     rows.append(row);save(rows)
     print(case,str(distance),direction,degree,cap,dec,row['stop_reason'],flush=True)
summary={'repository_commit':'1948ea20ae6d88395dc417c97581a78a427484da','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'reference_sha256':hashlib.sha256((BASE/'published_references.json').read_bytes()).hexdigest(),'python':platform.python_version(),'requests':len(rows),'decision_compatible':sum(r['decision_compatible'] for r in rows),'reference_contained':sum(r['reference_contained'] for r in rows),'errors':[r for r in rows if not r['decision_compatible'] or not r['reference_contained'] or r['returncode'] not in (0,2)],'grouped':[]}
for case in cases:
 for dist in ('1/100','1/1000','1/10000'):
  for cap in (64,256):
   for degree in (2,3):
    subset=[r for r in rows if r['case']==case and r['relative_distance']==dist and r['cap']==cap and r['degree']==degree]
    summary['grouped'].append({'case':case,'relative_distance':dist,'cap':cap,'degree':degree,'determined':sum(r['goal_met'] for r in subset),'total':len(subset)})
(OUT/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k not in ('grouped','errors')}), 'errors',len(summary['errors']),flush=True)
