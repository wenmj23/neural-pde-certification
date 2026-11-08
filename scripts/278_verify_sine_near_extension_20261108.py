"""Independent readback of post-hoc extended-mesh certificates."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results';D=B/'2026-11-08/sine_near_threshold_extension';OUT=B/'2026-11-08/sine_near_threshold_extension_readback.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(l,u,t):return 'within' if u<=t else 'exceeds' if l>t else 'unknown'
def main():
 if OUT.exists():raise FileExistsError(OUT)
 s=json.loads((D/'summary.json').read_text());rows=json.loads((D/'records.json').read_text());c=json.loads((ROOT/'configs/midterm_sine_decisions_2026-11-07.json').read_text());fixture=ROOT/c['fixture'];ref=ROOT/c['reference'];old=(ROOT/c['certificate_summary']).parent/'records.json'
 if sha(D/'records.json')!=s['records_sha256'] or sha(fixture)!=s['fixture_sha256'] or sha(ref)!=s['reference_sha256'] or sha(old)!=s['old_records_sha256']:raise ValueError('hash')
 for p,h in s['source_hashes_before_run'].items():
  if sha(ROOT/p)!=h:raise ValueError('source hash')
 nets={n['label']:n for n in json.loads(fixture.read_text())['networks']};refs={r['network']:r for r in json.loads(ref.read_text())};orig={(r['network'],r['degree'],r['intervals']):r for r in json.loads(old.read_text())};expected=set(product(nets,(2,3),(128,256)));seen=set();counts=Counter();matched=0
 for r in rows:
  k=(r['network'],r['degree'],r['mesh'])
  if k not in expected or k in seen or r['status']!='complete':raise ValueError('manifest')
  seen.add(k);q=refs[k[0]];rl,ru=F(q['lower']),F(q['upper']);l,u=F(r['lower']),F(r['upper']);prior=orig[k];same=(l,u)==(F(prior['lower']),F(prior['upper']))
  if not q['converged'] or r['frequency']!=nets[k[0]]['frequency'] or (F(r['reference_lower']),F(r['reference_upper']))!=(rl,ru) or not 0<=l<=rl<=ru<=u or not r['reference_contained'] or r['matches_old']!=same:raise ValueError('range')
  matched+=same
  for fac in (F(99,100),F(101,100)):
   tol=(rl+ru)*fac/2;v=dec(l,u,tol);w=dec(rl,ru,tol)
   if r['decisions'][str(fac)]!=v or (v!='unknown' and w!='unknown' and v!=w):raise ValueError('decision')
   counts[(k[1],k[2],str(fac),v)]+=1
 if seen!=expected or len(nets)!=12 or s['records']!=48 or s['complete']!=48 or s['contained']!=48 or s['matches_old']!=matched:raise ValueError('summary')
 out={'records':len(rows),'reference_contained':len(rows),'matches_old':matched,'decision_checks':96,'counts':[{'degree':k[0],'mesh':k[1],'factor':k[2],'decision':k[3],'count':v} for k,v in sorted(counts.items())],'records_sha256':sha(D/'records.json'),'scope':'post-exit data and saved-reference arithmetic; post-hoc thresholds and meshes, no timing or independent proof'}
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
