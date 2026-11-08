"""Post-exit sine decision readback and comparison to earlier certificates."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'results/2026-11-07/midterm_sine_decisions';OUT=ROOT/'results/2026-11-07/midterm_sine_decisions_readback.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(l,u,t):return 'within' if u<=t else 'exceeds' if l>t else 'unknown'
def main():
 if OUT.exists():raise FileExistsError(OUT)
 c=json.loads((D/'locked_config.json').read_text());s=json.loads((D/'summary.json').read_text());rows=json.loads((D/'records.json').read_text());fixture=ROOT/c['fixture'];ref=ROOT/c['reference'];old_dir=(ROOT/c['certificate_summary']).parent;old_summary=json.loads((old_dir/'summary.json').read_text())
 if sha(D/'records.json')!=s['records_sha256'] or sha(D/'locked_config.json')!=sha(ROOT/'configs/midterm_sine_decisions_2026-11-07.json') or sha(fixture)!=c['fixture_sha256'] or sha(ref)!=s['reference_sha256'] or sha(ref)!=old_summary['references_sha256'] or sha(old_dir/'records.json')!=old_summary['records_sha256']:raise ValueError('digest')
 for p,h in s['source_hashes_before_run'].items():
  if sha(ROOT/p)!=h:raise ValueError('source changed '+p)
 nets={n['label']:n for n in json.loads(fixture.read_text())['networks']};refs={q['network']:q for q in json.loads(ref.read_text())};old={(r['network'],r['degree'],r['intervals']):r for r in json.loads((old_dir/'records.json').read_text())};expected=set(product(nets,c['degrees'],c['meshes']));seen=set();stats=Counter();refstats=Counter();indexed={};matched=0
 for r in rows:
  key=(r['network'],r['degree'],r['mesh'])
  if key not in expected or key in seen or r['status']!='complete':raise ValueError('manifest')
  seen.add(key);n=nets[key[0]];q=refs[key[0]];rl,ru=F(q['lower']),F(q['upper']);L,U=F(r['lower']),F(r['upper'])
  if (r['frequency'],r['depth'],r['seed'])!=(n['frequency'],n['depth'],n['seed']) or not q['converged'] or F(r['reference_lower'])!=rl or F(r['reference_upper'])!=ru or not 0<=L<=rl<=ru<=U or not r['reference_contained']:raise ValueError('enclosure')
  matched+=F(old[key]['lower'])==L and F(old[key]['upper'])==U
  for t in c['tolerances']:
   v=dec(L,U,F(t));w=dec(rl,ru,F(t))
   if r['decisions'][t]!=v or r['reference_decisions'][t]!=w or (v!='unknown' and w!='unknown' and v!=w):raise ValueError('decision')
   stats[(r['degree'],r['mesh'],t,v)]+=1
  indexed[key]=r
 if seen!=expected or len(nets)!=12 or s['records']!=72 or s['complete']!=72 or s['reference_contained']!=72:raise ValueError('incomplete')
 for q in refs.values():
  for t in c['tolerances']:refstats[(t,dec(F(q['lower']),F(q['upper']),F(t)))]+=1
 transitions=Counter()
 for name,d,t in product(nets,c['degrees'],c['tolerances']):
  for a,b in ((16,32),(32,64)):transitions[(a,b,indexed[(name,d,a)]['decisions'][t],indexed[(name,d,b)]['decisions'][t])]+=1
 out={'records':72,'networks':12,'reference_contained':72,'decision_checks':288,'contradictions':0,'exactly_matched_old_certificates':matched,'reference_counts':[{'tolerance':k[0],'decision':k[1],'count':v} for k,v in sorted(refstats.items())],'groups':[{'degree':k[0],'mesh':k[1],'tolerance':k[2],'decision':k[3],'count':v} for k,v in sorted(stats.items())],'transitions':[{'from_mesh':k[0],'to_mesh':k[1],'from_decision':k[2],'to_decision':k[3],'count':v} for k,v in sorted(transitions.items())],'records_sha256':sha(D/'records.json'),'scope':'post-exit stored-value audit against old continuous numerical references; not independent proof'}
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps({k:v for k,v in out.items() if k!='groups'},indent=2))
if __name__=='__main__':main()
