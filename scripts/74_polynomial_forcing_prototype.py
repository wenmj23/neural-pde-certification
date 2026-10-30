"""Research prototype: polynomial source f, not public CLI or retraining.
- e''=r, exact zero boundaries; four prescribed forcing functions.
Independent known-solution mpmath interval reference is for checking only.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
from collections import Counter
import importlib.util,json,hashlib
from flint import arb
from arb_cubic_reconstruction import SixthOrderDeep,CubicVerifier
from arb_shaped_remainder import ShapedVerifier
from arb_certification import finite,exact_fraction
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)
spec=importlib.util.spec_from_file_location('mp_reference_24',ROOT/'scripts/24_depth_continuous_reference.py');mpmod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mpmod)
SOURCES={'constant_2':[2],'linear_6x':[0,6],'quadratic_12x2':[0,0,12],'quartic_30x4':[0,0,0,0,30]}
class PolynomialForcing(SixthOrderDeep):
 def __init__(self,path,label,sha,coeff,bits=103):super().__init__(path,label,sha,bits);self.coeff=coeff
 def source_derivative(self,x,j):
  out=arb(0)
  for i in reversed(range(j,len(self.coeff))):out=out*x+self.coeff[i]*factorial(i)//factorial(i-j)
  return finite(out,'polynomial source derivative')
 def residual(self,x):
  n=self.network_jet(x,2);return finite(-2*(self.q(x)*n[2]+(1-2*x)*n[1]-n[0])-self.source_derivative(x,0),'poly residual')
 def r_second(self,x):
  n=self.network_jet(x,4);return finite(-24*(self.q(x)*n[4]+(1-2*x)*n[3]-n[2])-self.source_derivative(x,2),'poly residual second')
 def r_third(self,x):
  n=self.network_jet(x,5);return finite(-120*(self.q(x)*n[5]+(1-2*x)*n[4]-n[3])-self.source_derivative(x,3),'poly residual third')
 def r_fourth(self,x):
  n=self.network_jet(x,6);return finite(-720*(self.q(x)*n[6]+(1-2*x)*n[5]-n[4])-self.source_derivative(x,4),'poly residual fourth')
class PolynomialReference(mpmod.MPDeepReference):
 def __init__(self,row,coeff):super().__init__(row);self.coeff=coeff
 def source(self,x):
  out=self.ctx.mpf(0)
  for c in reversed(self.coeff):out=out*x+c
  return out
 def error(self,x):
  if x not in self.node_cache:
   xx=self.point(x);n,_,_=self.network(xx);u=sum((c*(xx-xx**(i+2))/((i+1)*(i+2)) for i,c in enumerate(self.coeff)),self.ctx.mpf(0));self.node_cache[x]=xx*(1-xx)*n-u
  return self.node_cache[x]
 def residual(self,left,right):
  self.residual_calls+=1;c=self.ctx;x=c.mpf([self.point(left).a,self.point(right).b]);n,n1,n2=self.network(x,2)
  return -x*(1-x)*n2-2*(1-2*x)*n1+2*n-self.source(x)
def fr(d):return F(int(d['numerator']),int(d['denominator']))
fixture=ROOT/'configs/frozen_networks_holdout_2026-10-30.json';data=fixture.read_bytes();sha=hashlib.sha256(data).hexdigest();nets=json.loads(data)['networks'];rows=[];refs={}
for raw in nets:
 name=raw['label']
 for source,coeff in SOURCES.items():
  ref=PolynomialReference(raw,coeff).maximum(relative=F(1,10**7),absolute=F(1,10**13),cap=5000);refs[name+'__'+source]=ref;(OUT/'references.json').write_text(json.dumps(refs,indent=2)+'\n');rl,ru=fr(ref['lower']),fr(ref['upper']);print('REFERENCE',name,source,ref['converged'],ref['splits'],flush=True)
  for degree in (2,3):
   for n in (128,256):
    problem=PolynomialForcing(fixture,name,sha,coeff);v=ShapedVerifier(problem,'moment') if degree==2 else CubicVerifier(problem,'moment');v.uniform(n);r=v.certify();L,U=exact_fraction(r['lower']),exact_fraction(r['upper']);row={'network':name,'source':source,'source_coefficients':coeff,'degree':degree,'intervals':n,'lower':str(L),'upper':str(U),'remainder_upper':str(exact_fraction(r['remainder_upper'])),'reference_lower':str(rl),'reference_upper':str(ru),'reference_converged':ref['converged'],'reference_overlap':L<=ru and U>=rl,'reference_contained':L<=rl and U>=ru,'ratio':str(U/L) if L else None};rows.append(row);(OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n');print(name,source,degree,n,'contains',row['reference_contained'],flush=True)
summary={'fixture_sha256':sha,'networks':len(nets),'sources':len(SOURCES),'reference_cases':len(refs),'reference_converged':sum(x['converged'] for x in refs.values()),'certificates':len(rows),'reference_overlap':sum(x['reference_overlap'] for x in rows),'reference_contained':sum(x['reference_contained'] for x in rows),'gap_1_02':sum(F(x['ratio'])<=F(51,50) for x in rows if x['ratio']),'failures':[x for x in rows if not x['reference_overlap'] or not x['reference_converged']],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'existing sine-trained frozen networks transferred to polynomial PDE, no retraining; true solution for validation only'};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k!='failures'}),'failures',len(summary['failures']),flush=True)
