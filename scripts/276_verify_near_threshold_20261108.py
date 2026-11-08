"""Reconstruct post-hoc thresholds and decisions from original exact bounds."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results';D=B/'2026-11-08/near_threshold_posthoc';OUT=B/'2026-11-08/near_threshold_posthoc_readback.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(l,u,t):return 'within' if u<=t else 'exceeds' if l>t else 'unknown'
def main():
 if OUT.exists():raise FileExistsError(OUT)
 s=json.loads((D/'summary.json').read_text());rows=json.loads((D/'records.json').read_text())
 if sha(D/'records.json')!=s['records_sha256'] or sha(ROOT/'scripts/275_near_threshold_posthoc_20261108.py')!=s['script_sha256']:raise ValueError('hash')
 old={}
 for f in ('poly','sine'):
  p=B/f'2026-11-07/midterm_{f}_decisions/records.json'
  if sha(p)!=s['inputs_sha256'][f]:raise ValueError('input hash')
  for r in json.loads(p.read_text()):old[(f,r['network'],r['degree'],r['mesh'])]=r
 expected={(f,n,d,m,str(x)) for f,n,d,m in old for x in (F(99,100),F(101,100))};seen=set();groups=Counter();refs=Counter();conflicts=0
 for r in rows:
  key=(r['family'],r['network'],r['degree'],r['mesh'],r['factor'])
  if key not in expected or key in seen:raise ValueError('manifest')
  seen.add(key);p=old[key[:4]];l,u,rl,ru=map(F,(p['lower'],p['upper'],p['reference_lower'],p['reference_upper']));t=(rl+ru)*F(key[4])/2;v,w=dec(l,u,t),dec(rl,ru,t)
  if F(r['tolerance'])!=t or r['decision']!=v or r['reference_decision']!=w:raise ValueError('arithmetic')
  if v!='unknown' and w!='unknown' and v!=w:conflicts+=1
  groups[(key[0],key[2],key[3],key[4],v)]+=1;refs[(key[0],key[4],w)]+=1
 if seen!=expected or len(rows)!=432 or conflicts:raise ValueError('incomplete/conflict')
 out={'records':len(rows),'conflicts':conflicts,'groups':[{'family':k[0],'degree':k[1],'mesh':k[2],'factor':k[3],'decision':k[4],'count':v} for k,v in sorted(groups.items())],'reference_groups':[{'family':k[0],'factor':k[1],'reference_decision':k[2],'count':v} for k,v in sorted(refs.items())],'records_sha256':sha(D/'records.json'),'scope':'post-hoc reference-derived thresholds; exact stored arithmetic, not blind test or fresh certificates'}
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps({'records':len(rows),'conflicts':conflicts,'reference_groups':out['reference_groups']}))
if __name__=='__main__':main()
