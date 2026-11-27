import sys,json,itertools,pathlib,math
from fractions import Fraction as F
p=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(p/'scripts'));from flint import arb,ctx;from arb_reconstruction import norm_enclosure;from arb_certification import exact_fraction
ctx.prec=103;rows=[];errors=[]
def ball(q):return arb(q.numerator)/q.denominator
for n,a,b in itertools.product([4,5],[F(0),F(1,2),F(1)],[F(-1,8),F(0),F(1,8)]):
 power=[F(math.comb(n,k))*(-a)**(n-k)+(b if k==0 else 0) for k in range(n+1)];coeff=[sum((power[k]*F(math.comb(j,k),math.comb(n,k)) for k in range(j+1)),F(0)) for j in range(n+1)]
 norm=max(abs((x-a)**n+b) for x in {F(0),F(1),a});prev=None
 for depth in range(5):
  lo,hi=norm_enclosure(list(map(ball,coeff)),depth);l=exact_fraction(lo.lower());u=exact_fraction(hi.upper());ok=l<=norm<=u and (prev is None or l>=prev);prev=l
  rows.append({'degree':n,'a':str(a),'b':str(b),'depth':depth,'lower':str(l),'true_norm':str(norm),'upper':str(u),'valid':ok})
  if not ok:errors.append([n,str(a),str(b),depth])
controls={'depth_invalid':False}
try:norm_enclosure([arb(0),arb(1)],-1)
except ValueError:controls['depth_invalid']=True
out={'checks':len(rows),'valid':sum(x['valid'] for x in rows),'errors':errors,'controls':controls,'precision_bits':ctx.prec}
(p/'rows.json').write_text(json.dumps(rows,indent=2)+'\n');(p/'summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2));assert len(rows)==90 and not errors and all(controls.values())
