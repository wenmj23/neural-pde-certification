"""Independently recompute selected result-file counts and simple bound invariants.

Run with the path to the run_20261101 directory. Does not establish mathematical
correctness of the verifier or explain an earlier incomplete replay.
"""
import json, sys
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
root=Path(sys.argv[1]); errors=[]; checked={}
def check(ok,msg):
 if not ok: errors.append(msg)
def load(name):
 p=root/name/'records.json'; rows=json.loads(p.read_text()); s=json.loads((p.parent/'summary.json').read_text()); return rows,s
expected={'poly_fixed_tolerance':72,'matched_residual_ablation':192,'shaped_remainder_ablation':240,'matched_acceptance_benchmark':288,'holdout_benchmark':288,'analytic_shape_stress':120,'bump_mesh_comparison':36,'sine_mesh_comparison':72,'shifted_peak_comparison':48}
for name,n in expected.items():
 rows,s=load(name); check(len(rows)==n,f'{name}: expected {n}, got {len(rows)}');checked[name]=len(rows)
 for i,r in enumerate(rows):
  def c(ok,msg):check(ok,f'{name} row {i}: {msg}')
  c(r.get('error') is None,'exception')
  if name=='poly_fixed_tolerance':
   t=F(r['tolerance']);c(bool(r['trace']),'empty trace')
   for a in r['trace']:
    l,u=F(a['lower']),F(a['upper']);c(0<=l<=u,'invalid bounds');c(a['decision']==('within' if u<=t else 'exceeds' if l>t else 'unknown'),'decision differs')
   if r['trace']:c(r['decision']==r['trace'][-1]['decision'],'final differs')
  elif name in ('matched_residual_ablation','shaped_remainder_ablation','analytic_shape_stress'):
   if 'lower' not in r: c(False,'missing bounds');continue
   l,u=F(r['lower']),F(r['upper']);c(0<=l<=u,'invalid bounds')
   if name=='analytic_shape_stress':c(l<=F(r['true_maximum'])<=u,'exact max uncovered')
   if name=='matched_residual_ablation':c(l<=F(r['reference_upper']) and u>=F(r['reference_lower']),'reference disjoint')
  elif name in ('matched_acceptance_benchmark','holdout_benchmark'):
   t=F(r['tolerance']);c(bool(r['trace']),'empty trace')
   for a in r['trace']:c(a['within']==(F(a['upper'])<=t),'incorrect acceptance')
   c((r['status']=='within')==bool(r['trace'] and r['trace'][-1]['within']),'final status inconsistent')
  else:
   c(bool(r['trace']),'empty trace')
   for a in r['trace']:
    l,u=F(a['lower']),F(a['upper']);c(0<=l<=u,'invalid bounds');c(a['goal_met']==(l>0 and u<=F(51,50)*l),'invalid gap decision')
    if name in ('bump_mesh_comparison','sine_mesh_comparison'):c(l<=F(1,10000)<=u,'exact max uncovered')
 if name in ('matched_acceptance_benchmark','holdout_benchmark'):
  keys=[(r['network'],r['tolerance'],r['repeat'],r['method']) for r in rows];check(len(set(keys))==len(keys),f'{name}: duplicate request')
 if name=='poly_fixed_tolerance':
  keys=[(r['network'],r['degree'],r['tolerance']) for r in rows];check(len(set(keys))==len(keys),'poly: duplicate request');check(s['requests']==len(rows),'poly summary count mismatch');check(s['outcomes']==dict(Counter(r['decision'] for r in rows)),'poly outcome counts mismatch')
result={'checked_records':checked,'total_records':sum(checked.values()),'errors':errors,'ok':not errors,'scope':'counts, selected uniqueness and simple invariant checks; not a formal audit or proof'}
print(json.dumps(result,indent=2));sys.exit(0 if not errors else 1)
