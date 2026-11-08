"""Post-exit timing manifest, saved arithmetic and paired per-network summary."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from statistics import median
import hashlib,json,math
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'results/2026-11-08/sine_near_threshold_timing';OUT=ROOT/'results/2026-11-08/sine_near_threshold_timing_readback.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(l,u,t):return 'within' if u<=t else 'exceeds' if l>t else 'unknown'
def main():
 if OUT.exists():raise FileExistsError(OUT)
 c=json.loads((D/'locked_config.json').read_text());s=json.loads((D/'summary.json').read_text());rows=json.loads((D/'records.json').read_text());sc=json.loads((ROOT/c['source_config']).read_text());fixture=ROOT/sc['fixture'];ref=ROOT/sc['reference'];old=(ROOT/sc['certificate_summary']).parent/'records.json'
 if sha(D/'records.json')!=s['records_sha256'] or sha(D/'locked_config.json')!=sha(ROOT/'configs/sine_near_threshold_timing_2026-11-08.json'):raise ValueError('output hash')
 if any(sha(ROOT/p)!=h for p,h in s['source_hashes_before_run'].items()):raise ValueError('source hash')
 if {'fixture':sha(fixture),'reference':sha(ref),'old_records':sha(old)}!=s['inputs_sha256']:raise ValueError('input hash')
 nets={n['label']:n for n in json.loads(fixture.read_text())['networks']};refs={r['network']:r for r in json.loads(ref.read_text())};previous={(r['network'],r['degree'],r['intervals']):r for r in json.loads(old.read_text())};expected={(n,d,m,phase,i) for n,d,m in product(nets,c['degrees'],c['meshes']) for phase,count in (('warmup',c['warmups']),('timed',c['timed_repeats'])) for i in range(count)};seen=set();idx={}
 for r in rows:
  k=r['network'],r['degree'],r['mesh'],r['phase'],r['iteration']
  if k not in expected or k in seen:raise ValueError('manifest')
  seen.add(k);q=refs[k[0]];rl,ru=F(q['lower']),F(q['upper']);L,U=F(r['lower']),F(r['upper']);oldrow=previous[k[:3]];same=(L,U)==(F(oldrow['lower']),F(oldrow['upper']))
  if not q['converged'] or (F(r['reference_lower']),F(r['reference_upper']))!=(rl,ru) or not 0<=L<=rl<=ru<=U or not r['reference_contained'] or not same or r['matches_old']!=same or not math.isfinite(r['seconds']) or r['seconds']<=0:raise ValueError('range or timing')
  for factor in c['threshold_factors']:
   tol=(rl+ru)*F(factor)/2;v=dec(L,U,tol);w=dec(rl,ru,tol)
   if r['decisions'][factor]!=v or (v!='unknown' and w!='unknown' and v!=w):raise ValueError('decision')
  idx[k]=r
 if len(nets)!=12 or seen!=expected or len(rows)!=288 or s['records']!=288 or s['contained']!=288 or s['matches_old']!=288:raise ValueError('incomplete')
 grouped=[];ratios=[]
 for mesh,degree in product(c['meshes'],c['degrees']):
  medians=[median(idx[(name,degree,mesh,'timed',i)]['seconds'] for i in range(c['timed_repeats'])) for name in nets];counts=Counter(idx[(name,degree,mesh,'timed',0)]['decisions'][fac] for name in nets for fac in c['threshold_factors']);grouped.append({'mesh':mesh,'degree':degree,'median_of_network_medians_seconds':median(medians),'decision_counts':dict(counts)})
 for mesh in c['meshes']:
  medians=[median(idx[(n,3,mesh,'timed',i)]['seconds']/idx[(n,2,mesh,'timed',i)]['seconds'] for i in range(c['timed_repeats'])) for n in nets];ratios.append({'mesh':mesh,'networks':len(medians),'median_of_network_median_paired_ratios':median(medians),'cubic_faster_networks':sum(x<1 for x in medians)})
 out={'records':len(rows),'warmup':sum(r['phase']=='warmup' for r in rows),'timed':sum(r['phase']=='timed' for r in rows),'contained':len(rows),'matches_old':len(rows),'groups':grouped,'paired_cubic_over_quadratic':ratios,'records_sha256':sha(D/'records.json'),'scope':'post-exit completeness, recorded arithmetic and timings; no independent runtime reproducibility or mathematical proof'}
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
