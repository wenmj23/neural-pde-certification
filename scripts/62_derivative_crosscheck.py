"""High-precision point checks against interval jets; numerical diagnostics, not proof."""
from fractions import Fraction as F
from pathlib import Path
from math import factorial
import hashlib,json,time
import mpmath as mp
from flint import arb
from arb_cubic_reconstruction import SixthOrderDeep
from arb_reconstruction import point
from arb_certification import exact_fraction
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261030/derivative_crosscheck';OUT.mkdir(parents=True,exist_ok=False)
CASES=['depth1_k1_seed7','depth1_k5_seed17','depth2_k1_seed17','depth2_k5_seed27','depth3_k1_seed27','depth3_k5_seed7','depth1_k2_seed49','depth3_k2_seed49','depth3_k5_seed49']
mp.mp.dps=90;rows=[]
def mprat(q):return mp.mpf(q.numerator)/q.denominator
def make_fn(network):
 layers=network['layers']; converted=[]
 for layer in layers:
  weights=[[mp.mpf(float(w)) for w in row] for row in layer['weight']]; bias=[mp.mpf(float(b)) for b in layer['bias']];converted.append((weights,bias,layer['activation']))
 def v(x):
  state=[x]
  for weights,bias,activation in converted:
   state=[b+sum((w*y for w,y in zip(row,state)),mp.mpf('0')) for row,b in zip(weights,bias)]
   if activation=='tanh':state=[mp.tanh(y) for y in state]
  return x*(1-x)*state[0]
 return v
for name in CASES:
 inp=json.loads((ROOT.parent/'run_20261030/frozen_portfolio'/(name+'.case.json')).read_text());path=Path(inp['network']['fixture']); raw=next(r for r in json.loads(path.read_text())['networks'] if r['label']==name);fn=make_fn(raw);k=raw['frequency'];omega=k*mp.pi
 p=SixthOrderDeep(path,name,inp['network']['sha256'],103)
 for center in (F(1,7),F(1,2),F(6,7)):
  a,b=center-F(1,1000),center+F(1,1000);ball=point(a).lower().union(point(b).upper())
  ivs={0:p.residual(ball),2:p.r_second(ball),3:p.r_third(ball),4:p.r_fourth(ball)}
  for x in (a,center,b):
   xx=mprat(x);coeff=mp.taylor(fn,xx,6)
   for j in (0,2,3,4):
    expected=-factorial(j+2)*coeff[j+2]-omega**(j+2)*mp.sin(omega*xx+j*mp.pi/2)
    lower=exact_fraction(ivs[j].lower());upper=exact_fraction(ivs[j].upper());passed=mprat(lower)<=expected<=mprat(upper)
    rows.append({'name':name,'center':str(center),'point':str(x),'order':j,'arb_lower':str(lower),'arb_upper':str(upper),'mp_value':mp.nstr(expected,80),'passed':passed})
  (OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n')
 print(name,'complete',flush=True)
summary={'networks':len(CASES),'intervals_per_network':3,'points_per_interval':3,'derivative_orders':[0,2,3,4],'comparisons':len(rows),'passed':sum(x['passed'] for x in rows),'failures':[x for x in rows if not x['passed']],'mp_decimal_digits':90,'arb_precision_bits':103,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k!='failures'}),'failures',len(summary['failures']),flush=True)
