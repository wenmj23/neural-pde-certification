"""Post-exit recomputation of offline paths and summary from timing records."""
from collections import Counter
from itertools import product
from pathlib import Path
from statistics import median
import hashlib,json,math
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results/2026-11-08';SRC=B/'sine_near_threshold_timing/records.json';DATA=B/'sine_offline_strategy_comparison.json';OUT=B/'sine_offline_strategy_readback.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 if OUT.exists():raise FileExistsError(OUT)
 s=json.loads(DATA.read_text());rows=s['rows'];source=json.loads(SRC.read_text());cfg=json.loads((B/'sine_near_threshold_timing/locked_config.json').read_text())
 if sha(SRC)!=s['input_records_sha256'] or sha(ROOT/'scripts/281_sine_offline_strategy_comparison_20261108.py')!=s['script_sha256']:raise ValueError('hash')
 idx={(r['network'],r['degree'],r['mesh'],r['iteration']):r for r in source if r['phase']=='timed'};names=sorted({r['network'] for r in source});paths={'start64':[64,128,256],'start128':[128,256],'direct256':[256]};expected=set(product(names,cfg['degrees'],cfg['threshold_factors'],paths));seen=set();lookup={};counts=Counter()
 for r in rows:
  key=r['network'],r['degree'],r['factor'],r['strategy']
  if key not in expected or key in seen:raise ValueError('manifest')
  seen.add(key);times=[];sequences=[];decisions=[]
  for i in range(cfg['timed_repeats']):
   elapsed=0.0;used=[];decision='unknown'
   for mesh in paths[key[3]]:
    z=idx[(key[0],key[1],mesh,i)];elapsed+=z['seconds'];used.append(mesh);decision=z['decisions'][key[2]]
    if decision!='unknown':break
   times.append(elapsed);sequences.append(used);decisions.append(decision)
  if decisions!=[r['decision']]*cfg['timed_repeats'] or any(x!=r['meshes_used'] for x in sequences) or any(not math.isclose(a,b,rel_tol=1e-15,abs_tol=0) for a,b in zip(times,r['times'])) or not math.isclose(median(times),r['median_seconds'],rel_tol=1e-15):raise ValueError('arithmetic')
  counts[(key[1],key[3],r['decision'])]+=1;lookup[key]=r
 if seen!=expected or s['scenarios']!=144 or len(rows)!=144:raise ValueError('incomplete')
 stored={(x['degree'],x['strategy'],x['decision']):x['count'] for x in s['per_strategy_counts']}
 if stored!=counts:raise ValueError('count')
 for q in s['comparisons']:
  vals=[]
  for n,fac in product(names,cfg['threshold_factors']):
   x=lookup[(n,q['degree'],fac,q['numerator'])];y=lookup[(n,q['degree'],fac,q['denominator'])]
   if x['decision']!='unknown' and y['decision']!='unknown':
    if x['decision']!=y['decision']:raise ValueError('contradiction')
    vals.append(median(x['times'][i]/y['times'][i] for i in range(cfg['timed_repeats'])))
  if len(vals)!=q['matched_decisive_scenarios'] or sum(x<1 for x in vals)!=q['numerator_faster_scenarios'] or not math.isclose(median(vals),q['scenario_median_of_paired_ratios'],rel_tol=1e-15):raise ValueError('ratio')
 result={'records':len(seen),'counts_verified':True,'comparisons_verified':len(s['comparisons']),'data_sha256':sha(DATA),'scope':'post-exit recomputation of stored fixed-mesh timing arithmetic; not an executed online benchmark'}
 with OUT.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps(result))
if __name__=='__main__':main()
