"""Reusable pointwise polynomial-residual derivative diagnostic.

Supports original frozen-network sources and three fixed transfer sources.
Uses independent mpmath numerical differentiation at 60/90 decimal digits.
This diagnostic does not prove whole-cell bounds or certify floating inference.
"""
from pathlib import Path
from fractions import Fraction as F
from copy import deepcopy
import argparse,hashlib,importlib.util,json,os,sys
from mpmath.libmp import to_rational
from flint import arb
from arb_certification import exact_fraction
from arb_reconstruction import point
ROOT=Path(__file__).resolve().parents[1]
TRANSFER={'quadratic_12x2':[0,0,12],'cubic_20x3':[0,0,0,20],'quintic_42x5':[0,0,0,0,0,42]}
POINTS=(F(1,7),F(1,2),F(6,7));ORDERS=(0,2,3,4)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n').encode();tmp=p.with_suffix('.tmp')
 with tmp.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def point_reference():
 spec=importlib.util.spec_from_file_location('point_reference107',ROOT/'scripts/107_polynomial_derivative_crosscheck_20261101.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m.mp_value
def prototype():
 path=ROOT/'scripts/74_polynomial_forcing_prototype.py';source=path.read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
 if source.count(needle)!=1:raise ValueError('prototype layout changed')
 ns={'__file__':str(path)};exec(compile(source.replace(needle,'OUT=None'),'defs','exec'),ns);return ns['PolynomialForcing']
def run(fixtures,out,mode='all',omit_source_for_test=False):
 if mode not in ('original','transfer','all'):raise ValueError('unsupported source mode')
 out.mkdir(parents=True,exist_ok=False);point_value=point_reference();cls=prototype()
 if omit_source_for_test:
  class OmitSource(cls):
   def source_derivative(self,x,j):return arb(0)
  cls=OmitSource
 inputs=[];rows=[]
 for fixture in fixtures:
  fixture=fixture.resolve();digest=sha(fixture);nets=json.loads(fixture.read_bytes())['networks'];labels=[n['label'] for n in nets]
  if not labels or len(set(labels))!=len(labels):raise ValueError('fixture labels empty or not unique')
  inputs.append({'sha256':digest,'networks':len(nets)})
  for net in nets:
   sources={}
   if mode in ('original','all'):sources['original']=net['source_coefficients']
   if mode in ('transfer','all'):sources.update(TRANSFER)
   for source,coeff in sources.items():
    if not isinstance(coeff,list) or not coeff or any(type(c) is not int for c in coeff):raise ValueError('sources must use integer coefficients')
    p=cls(fixture,net['label'],digest,coeff); raw=deepcopy(net);raw['source_coefficients']=coeff;fns={0:p.residual,2:p.r_second,3:p.r_third,4:p.r_fourth}
    for x in POINTS:
     for order in ORDERS:
      ball=fns[order](point(x));lo,hi=exact_fraction(ball.lower()),exact_fraction(ball.upper());low=F(*to_rational(point_value(raw,x,order,60)._mpf_));high=F(*to_rational(point_value(raw,x,order,90)._mpf_));stable=abs(low-high)<=F(1,10**50)*max(F(1),abs(high));inside=lo<=high<=hi
      # At x>0, explicitly evaluate the exact forcing derivative; do not
      # infer it from a nominal polynomial degree, which could have zeros.
      from math import factorial
      forcing=sum((F(c)*factorial(i)//factorial(i-order)*x**(i-order) for i,c in enumerate(coeff) if i>=order),F(0));nonzero=forcing!=0
      expected_inside=not (omit_source_for_test and nonzero)
      rows.append({'fixture_sha256':digest,'network':net['label'],'source':source,'source_coefficients':coeff,'x':str(x),'order':order,'arb_lower':str(lo),'arb_upper':str(hi),'mp60':str(low),'mp90':str(high),'forcing_derivative':str(forcing),'stable':stable,'contained':inside,'expected_contained':expected_inside,'passed':stable and inside==expected_inside})
   print(net['label'],'completed',flush=True)
 atom(out/'records.json',rows);source_count=4 if mode=='all' else (3 if mode=='transfer' else 1);expected=sum(x['networks'] for x in inputs)*source_count*len(POINTS)*len(ORDERS)
 if len(rows)!=expected:raise ValueError('record count mismatch')
 s={'inputs':inputs,'mode':mode,'test_fault':'omit_source' if omit_source_for_test else 'none','transfer_sources':TRANSFER,'points':[str(x) for x in POINTS],'orders':ORDERS,'expected':expected,'records':len(rows),'stable':sum(r['stable'] for r in rows),'contained':sum(r['contained'] for r in rows),'expected_mismatch':sum(not r['expected_contained'] for r in rows),'passed':sum(r['passed'] for r in rows),'records_sha256':sha(out/'records.json'),'script_sha256':sha(Path(__file__)),'point_reference_sha256':sha(ROOT/'scripts/107_polynomial_derivative_crosscheck_20261101.py'),'failures':[{'fixture_sha256':r['fixture_sha256'],'network':r['network'],'source':r['source'],'x':r['x'],'order':r['order']} for r in rows if not r['passed']],'scope':'pointwise numerical diagnostic with exact binary64 weights; not an independent formal proof or continuous interval certificate'}
 atom(out/'summary.json',s);print(json.dumps(s),flush=True);return int(bool(s['failures']))
def main():
 p=argparse.ArgumentParser();p.add_argument('--fixtures',nargs='+',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--mode',choices=['original','transfer','all'],default='all');a=p.parse_args();return run(a.fixtures,a.out,a.mode)
if __name__=='__main__':sys.exit(main())
