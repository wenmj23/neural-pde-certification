"""Rational-polynomial test problems, separately generated rational Bernstein reference.
Synthetic residuals, not neural networks; no software-audit claim.
"""
from fractions import Fraction as F
from math import factorial,comb
from pathlib import Path
import json,hashlib,time
from flint import arb
from arb_problem_families import SmoothProblem
from arb_reconstruction import Verifier
from arb_shaped_remainder import ShapedVerifier
from arb_cubic_reconstruction import CubicVerifier
from arb_certification import exact_fraction
OUT=Path(__file__).resolve().parents[2]/'run_20261030/polynomial_audit';OUT.mkdir(parents=True,exist_ok=False)
# Coefficients of residual r in increasing powers; common scale 1/10000.
POLYS={'zero':[0],'constant':[1],'negative_constant':[-1],'linear':[0,1], 'quadratic':[1,-3,2],'cubic':[1,-7,15,-9],'node_hidden_quartic':[0,-2,11,-18,9], 'alternating_sextic':[1,-12,50,-100,100,-48,9], 'asymmetric_sextic':[-2,14,-33,27,0,-8,3]}
SCALE=F(1,10000)
def error_power(coeff):
 out=[F(0)]*(len(coeff)+2)
 for i,c in enumerate(coeff):
  v=F(c)*SCALE/F((i+1)*(i+2));out[1]+=v;out[i+2]-=v
 return out
# Entirely separate exact rational Bernstein implementation for reference.
def reference(coeff,depth=12):
 p=error_power(coeff); n=len(p)-1
 bern=[sum((p[j]*F(comb(k,j),comb(n,j)) for j in range(k+1)),F(0)) for k in range(n+1)]
 pieces=[bern]
 for _ in range(depth):
  nxt=[]
  for row in pieces:
   row=list(row);left=[row[0]];right=[row[-1]]
   while len(row)>1:
    row=[(a+b)/2 for a,b in zip(row[:-1],row[1:])];left.append(row[0]);right.append(row[-1])
   nxt.extend([left,list(reversed(right))])
  pieces=nxt
 lo=max(max(abs(p[0]),abs(p[-1])) for p in pieces); hi=max(abs(c) for p in pieces for c in p)
 return lo,hi
class Polynomial(SmoothProblem):
 def __init__(self,coeff):self.setup(103);self.coeff=coeff;self.name='rational_polynomial';self.boundary_error=(0,0)
 def derivative(self,x,n):
  self.guard();out=arb(0)
  for i in reversed(range(n,len(self.coeff))):out=out*x+arb(self.coeff[i]*factorial(i)//factorial(i-n))/10000
  return out
 def residual(self,x):return self.derivative(x,0)
 def r_second(self,x):return self.derivative(x,2)
 def r_third(self,x):return self.derivative(x,3)
 def r_fourth(self,x):return self.derivative(x,4)
meshes={'one':[F(0),F(1)],'four':[F(i,4) for i in range(5)],'sixteen':[F(i,16) for i in range(17)],'nonuniform':[F(0),F(1,17),F(1,3),F(7,10),F(1)]}
methods=[('quadratic_constant',lambda p:ShapedVerifier(p,'constant')),('quadratic_full',lambda p:ShapedVerifier(p,'full')),('quadratic_moment',lambda p:ShapedVerifier(p,'moment')),('cubic_constant',lambda p:CubicVerifier(p,'constant')),('cubic_moment',lambda p:CubicVerifier(p,'moment'))]
rows=[]; refs={}
for name,coeff in POLYS.items():
 lo,hi=reference(coeff);refs[name]={'lower':str(lo),'upper':str(hi),'depth':12,'known_exact':('0' if name=='zero' else str(SCALE/8) if name in ('constant','negative_constant') else None)}
 if name=='zero':assert lo==hi==0
 if name in ('constant','negative_constant'):assert lo==hi==SCALE/8
 for mesh,nodes in meshes.items():
  for method,construct in methods:
   t=time.perf_counter();problem=Polynomial(coeff);v=construct(problem);e=v.certify(nodes);lower=exact_fraction(e['lower']);upper=exact_fraction(e['upper']); compatible=lower<=hi and upper>=lo
   exact=refs[name]['known_exact'];exactok=(lower<=F(exact)<=upper) if exact is not None else None
   row={'case':name,'mesh':mesh,'method':method,'intervals':len(nodes)-1,'reference_lower':str(lo),'reference_upper':str(hi),'lower':str(lower),'upper':str(upper),'remainder_upper':str(exact_fraction(e['remainder_upper'])),'compatible':compatible,'exact_known_contained':exactok,'seconds':time.perf_counter()-t}
   rows.append(row);(OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n')
 print(name,'done',flush=True)
summary={'families':len(POLYS),'configurations':len(rows),'compatible':sum(r['compatible'] for r in rows),'exact_checks':sum(r['exact_known_contained'] is not None for r in rows),'exact_pass':sum(r['exact_known_contained'] is True for r in rows),'failures':[r for r in rows if not r['compatible'] or r['exact_known_contained'] is False],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(OUT/'reference.json').write_text(json.dumps(refs,indent=2)+'\n');(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k!='failures'}),'failures',len(summary['failures']),flush=True)
