"""Post-exit check repeated online paths, provenance and paired timings."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from statistics import median
import hashlib,json,math
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results/2026-11-08';D=B/'sine_fixed_tolerance_repeat';OUT=B/'sine_fixed_tolerance_repeat_readback.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(l,u,t):return 'within' if u<=t else 'exceeds' if l>t else 'unknown'
def main():
 if OUT.exists():raise FileExistsError(OUT)
 c=json.loads((D/'locked_config.json').read_text());s=json.loads((D/'summary.json').read_text());rows=json.loads((D/'records.json').read_text());bc=json.loads((ROOT/c['base_config']).read_text());sc=json.loads((ROOT/bc['source_config']).read_text());fixture=ROOT/sc['fixture'];ref=ROOT/sc['reference'];old=(ROOT/sc['certificate_summary']).parent/'records.json';original=B/'sine_fixed_tolerance_online/records.json'
 if sha(D/'records.json')!=s['records_sha256'] or sha(D/'locked_config.json')!=sha(ROOT/'configs/sine_fixed_tolerance_repeat_2026-11-08.json') or {'fixture':sha(fixture),'reference':sha(ref),'old_records':sha(old),'online_records':sha(original)}!=s['inputs_sha256']:raise ValueError('hash')
 if any(sha(ROOT/p)!=h for p,h in s['source_hashes_before_run'].items()):raise ValueError('source hash')
 nets={n['label']:n for n in json.loads(fixture.read_text())['networks']};refs={q['network']:q for q in json.loads(ref.read_text())};prev={(r['network'],r['degree'],r['intervals']):r for r in json.loads(old.read_text())};orig={(r['network'],r['degree'],r['tolerance'],r['strategy']):r for r in json.loads(original.read_text())};expected={(n,d,st,ph,i) for n,d,st in product(nets,bc['degrees'],bc['strategies']) for ph,num in (('warmup',c['warmups']),('timed',c['timed_repeats'])) for i in range(num)};seen=set();idx={};attempts=matches=paths=0;tol=F(c['tolerance'])
 for r in rows:
  key=r['network'],r['degree'],r['strategy'],r['phase'],r['iteration']
  if key not in expected or key in seen or r['status']!='complete':raise ValueError('manifest')
  seen.add(key);q=refs[key[0]];rl,ru=F(q['lower']),F(q['upper']);arr=r['attempts'];grid=bc['strategies'][key[2]]
  if (F(r['reference_lower']),F(r['reference_upper']))!=(rl,ru) or F(r['tolerance'])!=tol or not q['converged'] or not arr or [a['mesh'] for a in arr]!=grid[:len(arr)]:raise ValueError('request')
  for i,a in enumerate(arr):
   L,U=F(a['lower']),F(a['upper']);oldrow=prev[(key[0],key[1],a['mesh'])];same=(L,U)==(F(oldrow['lower']),F(oldrow['upper']));v=dec(L,U,tol);w=dec(rl,ru,tol)
   if not 0<=L<=rl<=ru<=U or not a['reference_contained'] or a['matches_old']!=same or a['decision']!=v or (v!='unknown' and w!='unknown' and v!=w) or (i<len(arr)-1 and v!='unknown') or not math.isfinite(a['seconds']) or a['seconds']<=0:raise ValueError('attempt')
   matches+=same;attempts+=1
  prior=orig[(key[0],key[1],str(tol),key[2])];samepath=arr[-1]['decision']==prior['decision'] and [a['mesh'] for a in arr]==[a['mesh'] for a in prior['attempts']]
  if not samepath or r['path_matches_previous']!=samepath or r['decision']!=arr[-1]['decision'] or (r['decision']=='unknown' and len(arr)!=len(grid)) or not math.isclose(r['core_seconds'],sum(a['seconds'] for a in arr),rel_tol=1e-12) or not math.isfinite(r['inclusive_seconds']) or r['inclusive_seconds']<r['core_seconds']:raise ValueError('stop/time')
  paths+=samepath;idx[key]=r
 if seen!=expected or len(nets)!=12 or len(rows)!=192 or s['records']!=192 or s['complete']!=192 or s['attempts']!=attempts or s['contained']!=attempts or s['matches_old']!=matches or s['paths_match_previous']!=paths:raise ValueError('summary')
 comparisons=[]
 for degree in bc['degrees']:
  network_ratios=[];faster_all=0
  for n in nets:
   ratios=[]
   for i in range(c['timed_repeats']):
    a=idx[(n,degree,'start16','timed',i)];b=idx[(n,degree,'start64','timed',i)]
    if a['decision']=='unknown' or b['decision']=='unknown' or a['decision']!=b['decision']:raise ValueError('missing/contradictory comparison')
    ratios.append(a['inclusive_seconds']/b['inclusive_seconds'])
   network_ratios.append(median(ratios));faster_all+=all(x<1 for x in ratios)
  comparisons.append({'degree':degree,'networks':len(network_ratios),'median_network_paired_time_ratio':median(network_ratios),'networks_median_faster_start16':sum(x<1 for x in network_ratios),'networks_all_three_faster_start16':faster_all})
 out={'records':len(rows),'warmup':sum(r['phase']=='warmup' for r in rows),'timed':sum(r['phase']=='timed' for r in rows),'attempts':attempts,'reference_contained':attempts,'matches_old':matches,'paths_match_previous':paths,'paired_start16_over_start64':comparisons,'records_sha256':sha(D/'records.json'),'scope':'post-exit arithmetic and recorded timing audit; one machine, repeated known networks; no formal proof'}
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
