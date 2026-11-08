"""Post-exit provenance, full manifest and exact decision arithmetic."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'results/2026-11-07/midterm_poly_decisions';OUT=ROOT/'results/2026-11-07/midterm_poly_decisions_readback.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def val(x):return F(int(x['numerator']),int(x['denominator']))
def dec(l,u,t):return 'within' if u<=t else 'exceeds' if l>t else 'unknown'
def main():
 if OUT.exists():raise FileExistsError(OUT)
 c=json.loads((D/'locked_config.json').read_text());s=json.loads((D/'summary.json').read_text());rows=json.loads((D/'records.json').read_text())
 if sha(D/'records.json')!=s['records_sha256'] or sha(D/'locked_config.json')!=sha(ROOT/'configs/midterm_poly_decisions_2026-11-07.json'):raise ValueError('digest')
 for p,h in s['source_hashes_before_run'].items():
  if sha(ROOT/p)!=h:raise ValueError('source changed '+p)
 cases={}
 for cfg in json.loads((ROOT/c['fixtures_config']).read_text())['fixtures']:
  fixture=ROOT/cfg['path'];reference=ROOT/cfg['reference']
  if sha(fixture)!=cfg['sha256'] or sha(reference)!=cfg['reference_sha256']:raise ValueError('fixture/reference hash')
  refs=json.loads(reference.read_text())
  for net in json.loads(fixture.read_text())['networks']:
   q=refs[net['label']];cases[net['label']]=(net,cfg['sha256'],val(q['lower']),val(q['upper']),q['converged'])
 expected=set(product(cases,c['degrees'],c['meshes']));seen=set();stats=Counter();refstats=Counter();indexed={};conflicts=[]
 for r in rows:
  key=(r['network'],r['degree'],r['mesh'])
  if key not in expected or key in seen or r['status']!='complete':raise ValueError('manifest')
  seen.add(key);net,h,rl,ru,converged=cases[key[0]];L,U=map(F,(r['lower'],r['upper']))
  if r['source']!=net['source_name'] or r['fixture_sha256']!=h or not converged or F(r['reference_lower'])!=rl or F(r['reference_upper'])!=ru or not 0<=L<=rl<=ru<=U or not r['reference_contained']:raise ValueError('enclosure')
  for t in c['tolerances']:
   decision=dec(L,U,F(t));reference=dec(rl,ru,F(t))
   if r['decisions'][t]!=decision or r['reference_decisions'][t]!=reference:raise ValueError('decision')
   stats[(r['degree'],r['mesh'],t,decision)]+=1
   if decision!='unknown' and reference!='unknown' and decision!=reference:conflicts.append((key,t))
  indexed[key]=r
 if seen!=expected or len(cases)!=24 or s['records']!=144 or s['complete']!=144 or s['reference_contained']!=144 or conflicts:raise ValueError('incomplete or contradictory')
 for name,(_,_,l,u,_) in cases.items():
  for t in c['tolerances']:refstats[(t,dec(l,u,F(t)))]+=1
 transitions=Counter()
 for name,d,t in product(cases,c['degrees'],c['tolerances']):
  for a,b in ((16,32),(32,64)):
   old=indexed[(name,d,a)]['decisions'][t];new=indexed[(name,d,b)]['decisions'][t];transitions[(a,b,old,new)]+=1
 out={'records':144,'networks':24,'reference_contained':144,'decision_checks':144*len(c['tolerances']),'contradictions':len(conflicts),'reference_counts':[{'tolerance':k[0],'decision':k[1],'count':v} for k,v in sorted(refstats.items())],'groups':[{'degree':k[0],'mesh':k[1],'tolerance':k[2],'decision':k[3],'count':v} for k,v in sorted(stats.items())],'transitions':[{'from_mesh':k[0],'to_mesh':k[1],'from_decision':k[2],'to_decision':k[3],'count':v} for k,v in sorted(transitions.items())],'records_sha256':sha(D/'records.json'),'scope':'post-exit record, input hash and exact-arithmetic audit; comparisons use saved continuous numerical references, not formal independent proof'}
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps({k:v for k,v in out.items() if k!='groups'},indent=2))
if __name__=='__main__':main()
