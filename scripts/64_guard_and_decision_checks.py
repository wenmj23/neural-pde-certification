"""Check fail-closed guards and exact boundary decisions."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
from flint import arb,ctx
import certification_pipeline as cp
from arb_reconstruction import Verifier
from arb_shaped_remainder import ShapedVerifier
from arb_cubic_reconstruction import CubicVerifier
from arb_certification import finite
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261030/guard_checks';OUT.mkdir(parents=True,exist_ok=False)
source=(ROOT/'scripts/61_exact_polynomial_audit.py').read_text().split('meshes=')[0];source=source.replace("OUT=Path(__file__).resolve().parents[2]/'run_20261030/polynomial_audit';OUT.mkdir(parents=True,exist_ok=False)","OUT=None");ns={'__file__':str(ROOT/'scripts/61_exact_polynomial_audit.py')};exec(compile(source,'polynomial_definitions','exec'),ns);P=ns['Polynomial']; rows=[]
def check(name,fn):
 try:ok=bool(fn());error=None
 except Exception as e:ok=False;error=f'{type(e).__name__}: {e}'
 rows.append({'name':name,'passed':ok,'error':error});print(name,ok,flush=True)
def raises(fn,cls=(ValueError,ArithmeticError,TypeError)):
 try:fn();return False
 except cls:return True
for typ in (ShapedVerifier,CubicVerifier):
 def create_bad(t=typ):
  p=P([1]);p.boundary_error=(1,0);t(p)
 check(typ.__name__+'_nonzero_boundary_rejected',lambda f=create_bad:raises(f))
 for label,nodes in [('missing_left',[F(1,10),F(1)]),('missing_right',[F(0),F(9,10)]),('repeated',[F(0),F(1,2),F(1,2),F(1)]),('reversed',[F(0),F(3,4),F(1,2),F(1)])]:
  check(typ.__name__+'_'+label+'_rejected',lambda t=typ,n=nodes:raises(lambda:t(P([1])).certify(n)))
 def changed_precision(t=typ):
  p=P([1]);v=t(p);ctx.prec=160
  try:return raises(lambda:v.certify())
  finally:ctx.prec=103
 check(typ.__name__+'_precision_change_rejected',changed_precision)
for label,body in [('duplicate_key','{"x":1,"x":2}'),('nan','{"x":NaN}'),('infinity','{"x":Infinity}')]:check('json_'+label,lambda b=body:raises(lambda:cp.strict_json(b)))
for invalid in [0,1.02,'nan','1/0','-1']:
 check('invalid_gap_'+str(invalid),lambda v=invalid:raises(lambda:cp.validate_settings({'gap_ratio':v})))
for invalid in [-1,0,79,701,103.0,True]:check('invalid_precision_'+str(invalid),lambda v=invalid:raises(lambda:cp.validate_settings({'precision_bits':v})))
s=cp.validate_settings({'goal':'decision','error_tolerance':'1/1000'})
for label,lo,hi,expected in [('upper_equal',F(0),F(1,1000),'within_tolerance'),('lower_equal',F(1,1000),F(2,1000),'undetermined'),('strict_above',F(1001,1000000),F(2,1000),'exceeds_tolerance'),('straddle',F(1,10000),F(2,1000),'undetermined')]:check('decision_'+label,lambda l=lo,h=hi,e=expected:cp.classify(l,h,s)['accuracy_decision']==e)
s0=cp.validate_settings({'goal':'decision','error_tolerance':'0'})
check('decision_exact_zero',lambda:cp.classify(F(0),F(0),s0)['accuracy_decision']=='within_tolerance')
check('decision_positive_vs_zero',lambda:cp.classify(F(1,1000),F(1,1000),s0)['accuracy_decision']=='exceeds_tolerance')
check('reject_reversed_enclosure',lambda:raises(lambda:cp.classify(F(2),F(1),s)))
check('reject_negative_enclosure',lambda:raises(lambda:cp.classify(F(-1),F(1),s)))
summary={'checks':len(rows),'passed':sum(x['passed'] for x in rows),'failed':[r for r in rows if not r['passed']],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n');(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',summary,flush=True)
