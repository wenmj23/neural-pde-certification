"""Audit saved fixed-tolerance sequential paths against original inputs."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from statistics import median
import hashlib,json,math
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results/2026-11-08';D=B/'sine_fixed_tolerance_online';OUT=B/'sine_fixed_tolerance_online_readback.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(l,u,t):return 'within' if u<=t else 'exceeds' if l>t else 'unknown'
def main():
 if OUT.exists():raise FileExistsError(OUT)
 c=json.loads((D/'locked_config.json').read_text());s=json.loads((D/'summary.json').read_text());rows=json.loads((D/'records.json').read_text());sc=json.loads((ROOT/c['source_config']).read_text());fixture=ROOT/sc['fixture'];ref=ROOT/sc['reference'];old=(ROOT/sc['certificate_summary']).parent/'records.json'
 if sha(D/'records.json')!=s['records_sha256'] or sha(D/'locked_config.json')!=sha(ROOT/'configs/sine_fixed_tolerance_online_2026-11-08.json') or {'fixture':sha(fixture),'reference':sha(ref),'old_records':sha(old)}!=s['inputs_sha256']:raise ValueError('digest')
 if any(sha(ROOT/p)!=h for p,h in s['source_hashes_before_run'].items()):raise ValueError('source hash')
 nets={n['label']:n for n in json.loads(fixture.read_text())['networks']};refs={q['network']:q for q in json.loads(ref.read_text())};prior={(r['network'],r['degree'],r['intervals']):r for r in json.loads(old.read_text())};expected=set(product(nets,c['degrees'],c['tolerances'],c['strategies']));seen=set();idx={};groups=Counter();attempt_count=0;matches=0
 for r in rows:
  key=r['network'],r['degree'],r['tolerance'],r['strategy']
  if key not in expected or key in seen or r['status']!='complete':raise ValueError('manifest')
  seen.add(key);q=refs[key[0]];rl,ru=F(q['lower']),F(q['upper']);tol=F(key[2]);arr=r['attempts'];grid=c['strategies'][key[3]]
  if not q['converged'] or r['frequency']!=nets[key[0]]['frequency'] or (F(r['reference_lower']),F(r['reference_upper']))!=(rl,ru) or r['reference_decision']!=dec(rl,ru,tol) or not arr or [a['mesh'] for a in arr]!=grid[:len(arr)]:raise ValueError('request')
  for i,a in enumerate(arr):
   L,U=F(a['lower']),F(a['upper']);prev=prior[(key[0],key[1],a['mesh'])];same=(L,U)==(F(prev['lower']),F(prev['upper']));v=dec(L,U,tol);w=dec(rl,ru,tol)
   if not 0<=L<=rl<=ru<=U or not a['reference_contained'] or a['matches_old']!=same or a['decision']!=v or (v!='unknown' and w!='unknown' and v!=w) or (i<len(arr)-1 and v!='unknown') or not math.isfinite(a['seconds']) or a['seconds']<=0:raise ValueError('attempt')
   if i and (a['residual_calls_cumulative']<arr[i-1]['residual_calls_cumulative'] or a['derivative_calls_cumulative']<arr[i-1]['derivative_calls_cumulative']):raise ValueError('counters')
   matches+=same;attempt_count+=1
  if r['decision']!=arr[-1]['decision'] or (r['decision']=='unknown' and len(arr)!=len(grid)) or not math.isclose(r['core_seconds'],sum(a['seconds'] for a in arr),rel_tol=1e-12) or not math.isfinite(r['inclusive_seconds']) or r['inclusive_seconds']<r['core_seconds']:raise ValueError('stop/time')
  groups[(key[1],key[2],key[3],r['decision'])]+=1;idx[key]=r
 if seen!=expected or len(rows)!=96 or s['records']!=96 or s['complete']!=96 or s['attempts']!=attempt_count or s['contained']!=attempt_count or s['matches_old']!=matches:raise ValueError('summary')
 pair=[]
 for d,t in product(c['degrees'],c['tolerances']):
  vals=[];status=Counter()
  for name in nets:
   a,b=idx[(name,d,t,'start16')],idx[(name,d,t,'start64')];status[(a['decision'],b['decision'])]+=1
   if a['decision']!='unknown' and b['decision']!='unknown':
    if a['decision']!=b['decision']:raise ValueError('contradiction')
    vals.append(a['inclusive_seconds']/b['inclusive_seconds'])
  pair.append({'degree':d,'tolerance':t,'matched_decisive_networks':len(vals),'median_one_run_start16_over_start64':median(vals) if vals else None,'start16_faster':sum(x<1 for x in vals),'paired_outcomes':[{'start16':k[0],'start64':k[1],'count':v} for k,v in sorted(status.items())]})
 out={'records':len(rows),'networks':len(nets),'attempts':attempt_count,'reference_contained':attempt_count,'matches_old':matches,'groups':[{'degree':k[0],'tolerance':k[1],'strategy':k[2],'decision':k[3],'count':v} for k,v in sorted(groups.items())],'paired_one_run':pair,'records_sha256':sha(D/'records.json'),'scope':'post-exit saved arithmetic, provenance and stop checks; single-run timing and numerical reference, not mathematical proof'}
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
