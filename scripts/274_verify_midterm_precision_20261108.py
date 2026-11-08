"""Post-exit audit of stored precision-sensitivity records."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results';D=B/'2026-11-08/midterm_precision_check';OUT=B/'2026-11-08/midterm_precision_check_readback.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def decision(l,u,t):return 'within' if u<=t else 'exceeds' if l>t else 'unknown'
def main():
 if OUT.exists():raise FileExistsError(OUT)
 c=json.loads((D/'locked_config.json').read_text());s=json.loads((D/'summary.json').read_text());rows=json.loads((D/'records.json').read_text())
 if sha(D/'records.json')!=s['records_sha256'] or sha(D/'locked_config.json')!=sha(ROOT/'configs/midterm_precision_check_2026-11-08.json'):raise ValueError('digest')
 for p,h in s['source_hashes_before_run'].items():
  if sha(ROOT/p)!=h:raise ValueError('source '+p)
 old={}
 for fam in ('poly','sine'):
  path=B/f'2026-11-07/midterm_{fam}_decisions/records.json'
  if sha(path)!=s['previous_records_sha256'][fam]:raise ValueError('prior hash')
  old[fam]={(r['network'],r['degree'],r['mesh']):r for r in json.loads(path.read_text())}
 expected={(f,n,b,d,m) for f in ('poly','sine') for n,b,d,m in product(c[f+'_networks'],c['bits'],c['degrees'],c['meshes'])};seen=set();idx={};outcomes=Counter();matched=0
 for r in rows:
  key=r['family'],r['network'],r['bits'],r['degree'],r['mesh']
  if key not in expected or key in seen or r['status']!='complete':raise ValueError('manifest')
  seen.add(key);p=old[key[0]][(key[1],key[3],key[4])];L,U,rl,ru,tol=map(F,(r['lower'],r['upper'],r['reference_lower'],r['reference_upper'],r['tolerance']))
  if (rl,ru)!=(F(p['reference_lower']),F(p['reference_upper'])) or tol!=F(c[key[0]+'_tolerance']) or not 0<=L<=rl<=ru<=U or not r['reference_contained']:raise ValueError('bounds/reference')
  dec=decision(L,U,tol);ref=decision(rl,ru,tol)
  if r['decision']!=dec or r['reference_decision']!=ref or (dec!='unknown' and ref!='unknown' and dec!=ref):raise ValueError('decision')
  if key[2]==103:
   same=(L,U)==(F(p['lower']),F(p['upper']))
   if same!=r['matches_prior_103']:raise ValueError('prior match')
   matched+=same
  idx[key]=r;outcomes[(key[0],key[2],dec)]+=1
 if seen!=expected or len(rows)!=96 or s['complete']!=96 or s['contained']!=96 or matched!=s['matches_prior_103']:raise ValueError('incomplete')
 widths=Counter();changes=Counter();max_rel=F(0)
 for f in ('poly','sine'):
  for n,d,m in product(c[f+'_networks'],c['degrees'],c['meshes']):
   for a,b in ((80,103),(103,160)):
    x,y=idx[(f,n,a,d,m)],idx[(f,n,b,d,m)];l1,u1,l2,u2=map(F,(x['lower'],x['upper'],y['lower'],y['upper']));w1,w2=u1-l1,u2-l2
    relation='narrower' if w2<w1 else 'equal' if w2==w1 else 'wider';widths[(f,a,b,relation)]+=1;changes[(f,a,b,x['decision'],y['decision'])]+=1
    if w2>w1 and w1>0:max_rel=max(max_rel,(w2-w1)/w1)
 out={'records':len(rows),'reference_contained':len(rows),'matched_103_prior':matched,'width_comparisons':[{'family':k[0],'from_bits':k[1],'to_bits':k[2],'comparison':k[3],'count':v} for k,v in sorted(widths.items())],'decision_transitions':[{'family':k[0],'from_bits':k[1],'to_bits':k[2],'from_decision':k[3],'to_decision':k[4],'count':v} for k,v in sorted(changes.items())],'max_relative_width_increase':str(max_rel),'max_relative_width_increase_approx':float(max_rel),'records_sha256':sha(D/'records.json'),'scope':'post-exit stored rational arithmetic and prior reference comparison; not an independent interval-implementation proof'}
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
