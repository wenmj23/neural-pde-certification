"""Post-exit manifest, exact arithmetic and report-table audit."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'results/2026-11-07/analytic_report_controls';OUT=ROOT/'results/2026-11-07/analytic_report_controls_readback.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 if OUT.exists():raise FileExistsError(OUT)
 c=json.loads((D/'locked_config.json').read_text());s=json.loads((D/'summary.json').read_text());rows=json.loads((D/'records.json').read_text())
 if sha(D/'records.json')!=s['records_sha256'] or sha(D/'locked_config.json')!=s['config_sha256'] or sha(ROOT/'configs/analytic_report_controls_2026-11-07.json')!=s['config_sha256'] or sha(ROOT/'scripts/263_analytic_report_controls_20261107.py')!=s['script_sha256']:raise ValueError('input/output digest')
 expected={('sine',k,m,d) for k,m,d in product(c['frequencies'],c['meshes'],c['degrees'])}|{('quartic',None,m,d) for m,d in product(c['quartic_meshes'],c['degrees'])};seen=set();stats=Counter();idx={}
 for r in rows:
  key=(r['family'],r['frequency'],r['mesh'],r['degree']);
  if key not in expected or key in seen:raise ValueError('manifest')
  seen.add(key);L,U,rem,K=map(F,(r['lower'],r['upper'],r['remainder_upper'],r['known_error']));contain=L<=K<=U;ratio=L>0 and U<=F(c['ratio_target'])*L
  if K!=F(1,10000) or not 0<=L<=U or rem<0 or r['contains_exact']!=contain or not contain or r['ratio_met']!=ratio:raise ValueError('arithmetic')
  for t in c['tolerances']:
   T=F(t);decision='within' if U<=T else 'exceeds' if L>T else 'undetermined'
   if r['decisions'][t]!=decision or (decision=='within' and K>T) or (decision=='exceeds' and K<=T):raise ValueError('decision')
   stats[(r['degree'],t,decision)]+=1
  if r['family']=='quartic' and rem!=0:raise ValueError('expected zero residual interpolation remainder')
  idx[key]=r
 if seen!=expected or len(rows)!=s['records'] or s['expected']!=78 or s['contains_exact']!=78:raise ValueError('incomplete')
 mesh_table=[]
 for k in c['frequencies']:
  row={'frequency':k}
  for d in c['degrees']:
   row['degree_'+str(d)]=next((m for m in c['meshes'] if idx[('sine',k,m,d)]['ratio_met']),None)
  mesh_table.append(row)
 out={'records':len(rows),'contained':78,'quartic_zero_remainder':6,'mesh_first_ratio_goal':mesh_table,'decision_counts':[{'degree':k[0],'tolerance':k[1],'decision':k[2],'count':v} for k,v in sorted(stats.items())],'records_sha256':sha(D/'records.json'),'source_hashes_after_run':{p:sha(ROOT/p) for p in ('scripts/arb_problem_families.py','scripts/arb_cubic_reconstruction.py','scripts/arb_shaped_remainder.py','scripts/arb_reconstruction.py','scripts/arb_certification.py')},'scope':'independent post-exit manifest and exact stored-value arithmetic; not a separate implementation of the verifier'}
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
