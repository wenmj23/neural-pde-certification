"""Post-hoc arithmetic replay of fixed-mesh timings, not an online benchmark."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
from statistics import median
import hashlib,json
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results/2026-11-08';D=B/'sine_near_threshold_timing';OUT=B/'sine_offline_strategy_comparison.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 if OUT.exists():raise FileExistsError(OUT)
 audit=json.loads((B/'sine_near_threshold_timing_readback.json').read_text());p=D/'records.json'
 if sha(p)!=audit['records_sha256']:raise ValueError('source digest')
 c=json.loads((D/'locked_config.json').read_text());rows=json.loads(p.read_text());idx={(r['network'],r['degree'],r['mesh'],r['iteration']):r for r in rows if r['phase']=='timed'};names=sorted({r['network'] for r in rows});paths={'start64':[64,128,256],'start128':[128,256],'direct256':[256]};out=[]
 for n in names:
  for d in c['degrees']:
   for fac in c['threshold_factors']:
    for strategy,meshes in paths.items():
     outcomes=[];times=[];sequences=[]
     for i in range(c['timed_repeats']):
      elapsed=0.0;used=[];decision='unknown'
      for m in meshes:
       r=idx[(n,d,m,i)];elapsed+=r['seconds'];used.append(m);decision=r['decisions'][fac]
       if decision!='unknown':break
      times.append(elapsed);outcomes.append(decision);sequences.append(used)
     if len(set(outcomes))!=1 or any(x!=sequences[0] for x in sequences):raise ValueError('repeat decision mismatch')
     out.append({'network':n,'degree':d,'factor':fac,'strategy':strategy,'decision':outcomes[0],'meshes_used':sequences[0],'times':times,'median_seconds':median(times)})
 counts=Counter((r['degree'],r['strategy'],r['decision']) for r in out);comparisons=[]
 for d in c['degrees']:
  lookup={(r['network'],r['factor'],r['strategy']):r for r in out if r['degree']==d}
  for a,b in (('start64','start128'),('start128','direct256'),('start64','direct256')):
   eligible=[]
   for n in names:
    for fac in c['threshold_factors']:
     x,y=lookup[(n,fac,a)],lookup[(n,fac,b)]
     if x['decision']!='unknown' and y['decision']!='unknown':
      if x['decision']!=y['decision']:raise ValueError('inconsistent results')
      eligible.append(median(x['times'][i]/y['times'][i] for i in range(c['timed_repeats'])))
   comparisons.append({'degree':d,'numerator':a,'denominator':b,'matched_decisive_scenarios':len(eligible),'scenario_median_of_paired_ratios':median(eligible) if eligible else None,'numerator_faster_scenarios':sum(x<1 for x in eligible)})
 result={'scenarios':len(out),'per_strategy_counts':[{'degree':k[0],'strategy':k[1],'decision':k[2],'count':v} for k,v in sorted(counts.items())],'comparisons':comparisons,'input_records_sha256':sha(p),'script_sha256':sha(Path(__file__)),'rows':out,'scope':'post-hoc sum of separately measured cold-object certify() times; omits object construction, I/O and cache reuse; thresholds depend on known references; not executed online strategy or matched total-budget experiment'}
 with OUT.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
if __name__=='__main__':main()
