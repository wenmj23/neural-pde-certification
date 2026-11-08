"""Read back actual adaptive paths, stopping rule, bounds and timing."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from statistics import median
import hashlib,json,math
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results/2026-11-08';D=B/'sine_online_near_threshold';OUT=B/'sine_online_near_threshold_readback.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(l,u,t):return 'within' if u<=t else 'exceeds' if l>t else 'unknown'
def main():
 if OUT.exists():raise FileExistsError(OUT)
 c=json.loads((D/'locked_config.json').read_text());s=json.loads((D/'summary.json').read_text());rows=json.loads((D/'records.json').read_text());sc=json.loads((ROOT/c['source_config']).read_text());fixture=ROOT/sc['fixture'];ref=ROOT/sc['reference'];old=(ROOT/sc['certificate_summary']).parent/'records.json'
 if sha(D/'records.json')!=s['records_sha256'] or sha(D/'locked_config.json')!=sha(ROOT/'configs/sine_online_near_threshold_2026-11-08.json') or {'fixture':sha(fixture),'reference':sha(ref),'old_records':sha(old)}!=s['inputs_sha256']:raise ValueError('hash')
 if any(sha(ROOT/p)!=h for p,h in s['source_hashes_before_run'].items()):raise ValueError('source hash')
 nets={n['label']:n for n in json.loads(fixture.read_text())['networks']};refs={q['network']:q for q in json.loads(ref.read_text())};previous={(r['network'],r['degree'],r['intervals']):r for r in json.loads(old.read_text())};off=json.loads((B/'sine_offline_strategy_comparison.json').read_text());offaudit=json.loads((B/'sine_offline_strategy_readback.json').read_text())
 if sha(B/'sine_offline_strategy_comparison.json')!=offaudit['data_sha256']:raise ValueError('offline digest')
 offline={(r['network'],r['degree'],r['factor'],r['strategy']):r for r in off['rows']};expected=set(product(nets,c['degrees'],c['threshold_factors'],c['strategies']));seen=set();idx={};counts=Counter();attempt_count=0;matches=0;path_match=0
 for r in rows:
  key=r['network'],r['degree'],r['factor'],r['strategy']
  if key not in expected or key in seen or r['status']!='complete':raise ValueError('manifest')
  seen.add(key);q=refs[key[0]];rl,ru=F(q['lower']),F(q['upper']);tol=(rl+ru)*F(key[2])/2;arr=r['attempts'];meshes=c['strategies'][key[3]]
  if not q['converged'] or (F(r['reference_lower']),F(r['reference_upper']),F(r['tolerance']))!=(rl,ru,tol) or not arr or [a['mesh'] for a in arr]!=meshes[:len(arr)]:raise ValueError('request')
  for i,a in enumerate(arr):
   L,U=F(a['lower']),F(a['upper']);prev=previous[(key[0],key[1],a['mesh'])];same=(L,U)==(F(prev['lower']),F(prev['upper']));v=dec(L,U,tol);w=dec(rl,ru,tol)
   if not 0<=L<=rl<=ru<=U or not a['reference_contained'] or a['matches_old']!=same or a['decision']!=v or (v!='unknown' and w!='unknown' and v!=w) or (i<len(arr)-1 and v!='unknown') or not math.isfinite(a['seconds']) or a['seconds']<=0:raise ValueError('attempt')
   if i and (a['residual_calls_cumulative']<arr[i-1]['residual_calls_cumulative'] or a['derivative_calls_cumulative']<arr[i-1]['derivative_calls_cumulative']):raise ValueError('counters')
   matches+=same;attempt_count+=1
  if r['decision']!=arr[-1]['decision'] or (r['decision']=='unknown' and len(arr)!=len(meshes)) or not math.isclose(r['core_seconds'],sum(a['seconds'] for a in arr),rel_tol=1e-12) or not math.isfinite(r['inclusive_seconds']) or r['inclusive_seconds']<r['core_seconds']:raise ValueError('stop/time')
  prior=offline[key];path_match+=(r['decision']==prior['decision'] and [a['mesh'] for a in arr]==prior['meshes_used']);counts[(key[1],key[3],r['decision'])]+=1;idx[key]=r
 if seen!=expected or len(rows)!=96 or s['records']!=96 or s['complete']!=96 or s['attempts']!=attempt_count or s['contained']!=attempt_count or s['matches_old']!=matches:raise ValueError('summary')
 paired=[]
 for d in c['degrees']:
  vals=[]
  for n,fac in product(nets,c['threshold_factors']):
   a,b=idx[(n,d,fac,'start64')],idx[(n,d,fac,'start128')]
   if a['decision']!='unknown' and b['decision']!='unknown':
    if a['decision']!=b['decision']:raise ValueError('strategy conflict')
    vals.append(a['inclusive_seconds']/b['inclusive_seconds'])
  paired.append({'degree':d,'matched_decisive_scenarios':len(vals),'median_inclusive_start64_over_start128':median(vals) if vals else None,'start64_faster':sum(x<1 for x in vals)})
 out={'records':len(rows),'attempts':attempt_count,'contained':attempt_count,'matches_old':matches,'offline_paths_matching':path_match,'outcomes':[{'degree':k[0],'strategy':k[1],'decision':k[2],'count':v} for k,v in sorted(counts.items())],'one_run_paired_inclusive_times':paired,'records_sha256':sha(D/'records.json'),'scope':'post-exit recorded arithmetic and stopping checks; one randomized run, no new mathematical proof'}
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
