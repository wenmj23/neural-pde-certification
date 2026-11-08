"""Readback of polynomial early-stop paths against fixed input references."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from statistics import median
import hashlib,json,math
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results/2026-11-08';D=B/'poly_fixed_tolerance_online';OUT=B/'poly_fixed_tolerance_online_readback.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def val(x):return F(int(x['numerator']),int(x['denominator']))
def dec(l,u,t):return 'within' if u<=t else 'exceeds' if l>t else 'unknown'
def main():
 if OUT.exists():raise FileExistsError(OUT)
 c=json.loads((D/'locked_config.json').read_text());s=json.loads((D/'summary.json').read_text());rows=json.loads((D/'records.json').read_text());bc=json.loads((ROOT/c['source_config']).read_text());prior=ROOT/'results/2026-11-07/midterm_poly_decisions/records.json'
 if sha(D/'records.json')!=s['records_sha256'] or sha(D/'locked_config.json')!=sha(ROOT/'configs/poly_fixed_tolerance_online_2026-11-08.json') or sha(prior)!=s['prior_records_sha256']:raise ValueError('hash')
 if any(sha(ROOT/p)!=h for p,h in s['source_hashes_before_run'].items()):raise ValueError('source hash')
 old={(r['network'],r['degree'],r['mesh']):r for r in json.loads(prior.read_text())};nets={}
 for f in json.loads((ROOT/bc['fixtures_config']).read_text())['fixtures']:
  path,ref=ROOT/f['path'],ROOT/f['reference']
  if sha(path)!=f['sha256'] or sha(ref)!=f['reference_sha256']:raise ValueError('fixture/reference')
  refs=json.loads(ref.read_text())
  for n in json.loads(path.read_text())['networks']:
   q=refs[n['label']]
   if n['label'] in nets or not q['converged']:raise ValueError('network/ref')
   nets[n['label']]=(n,f['sha256'],val(q['lower']),val(q['upper']))
 expected=set(product(nets,c['degrees'],c['tolerances'],c['strategies']));seen=set();idx={};groups=Counter();attempts=matches=0;stopmeshes=Counter()
 for r in rows:
  key=r['network'],r['degree'],r['tolerance'],r['strategy']
  if key not in expected or key in seen or r['status']!='complete':raise ValueError('manifest')
  seen.add(key);n,h,rl,ru=nets[key[0]];tol=F(key[2]);arr=r['attempts'];grid=c['strategies'][key[3]]
  if (r['source'],r['fixture_sha256'])!=(n['source_name'],h) or (F(r['reference_lower']),F(r['reference_upper']))!=(rl,ru) or r['reference_decision']!=dec(rl,ru,tol) or not arr or [a['mesh'] for a in arr]!=grid[:len(arr)]:raise ValueError('request')
  for i,a in enumerate(arr):
   L,U=F(a['lower']),F(a['upper']);p=old[(key[0],key[1],a['mesh'])];same=(L,U)==(F(p['lower']),F(p['upper']));v,w=dec(L,U,tol),dec(rl,ru,tol)
   if not 0<=L<=rl<=ru<=U or not a['reference_contained'] or a['matches_old']!=same or a['decision']!=v or (v!='unknown' and w!='unknown' and v!=w) or (i<len(arr)-1 and v!='unknown') or not math.isfinite(a['seconds']) or a['seconds']<=0:raise ValueError('attempt')
   if i and (a['residual_calls_cumulative']<arr[i-1]['residual_calls_cumulative'] or a['derivative_calls_cumulative']<arr[i-1]['derivative_calls_cumulative']):raise ValueError('counters')
   matches+=same;attempts+=1
  if r['decision']!=arr[-1]['decision'] or (r['decision']=='unknown' and len(arr)!=len(grid)) or not math.isclose(r['core_seconds'],sum(a['seconds'] for a in arr),rel_tol=1e-12) or not math.isfinite(r['inclusive_seconds']) or r['inclusive_seconds']<r['core_seconds']:raise ValueError('stop/time')
  groups[(key[1],key[2],key[3],r['decision'])]+=1;stopmeshes[(key[1],key[2],key[3],arr[-1]['mesh'])]+=1;idx[key]=r
 if seen!=expected or len(nets)!=24 or len(rows)!=192 or s['records']!=192 or s['complete']!=192 or s['attempts']!=attempts or s['contained']!=attempts or s['matches_old']!=matches:raise ValueError('summary')
 paired=[]
 for d,t in product(c['degrees'],c['tolerances']):
  ratios=[];outcomes=Counter()
  for n in nets:
   a,b=idx[(n,d,t,'start16')],idx[(n,d,t,'start64')];outcomes[(a['decision'],b['decision'])]+=1
   if a['decision']!='unknown' and b['decision']!='unknown':
    if a['decision']!=b['decision']:raise ValueError('contradiction')
    ratios.append(a['inclusive_seconds']/b['inclusive_seconds'])
  paired.append({'degree':d,'tolerance':t,'matched_decisive_networks':len(ratios),'median_one_run_start16_over_start64':median(ratios) if ratios else None,'start16_faster':sum(v<1 for v in ratios),'paired_outcomes':[{'start16':k[0],'start64':k[1],'count':v} for k,v in sorted(outcomes.items())]})
 out={'records':len(rows),'networks':len(nets),'attempts':attempts,'reference_contained':attempts,'matches_old':matches,'groups':[{'degree':k[0],'tolerance':k[1],'strategy':k[2],'decision':k[3],'count':v} for k,v in sorted(groups.items())],'stop_mesh_counts':[{'degree':k[0],'tolerance':k[1],'strategy':k[2],'mesh':k[3],'count':v} for k,v in sorted(stopmeshes.items())],'paired_one_run':paired,'records_sha256':sha(D/'records.json'),'scope':'post-exit saved arithmetic, provenance, stop checks and single-run timing; not mathematical proof'}
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps({k:v for k,v in out.items() if k!='groups'},indent=2))
if __name__=='__main__':main()
