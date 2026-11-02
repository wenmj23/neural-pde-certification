"""Cross-check full-cell derivative enclosures with independently propagated interval jets.

Mpmath interval tanh is built from exp, and its derivative jets use y'=z'(1-y²).
This differs from Arb's tanh derivative polynomials, but shares the same input,
PDE and mathematical premises. Sample values are diagnostic, not a proof.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import argparse, hashlib, importlib.util, json, os, sys
import mpmath as mp
from flint import arb
from arb_certification import exact_fraction
from arb_reconstruction import point

ROOT=Path(__file__).resolve().parents[1]
IV=mp.iv
IV.dps=50
ORDERS=(0,2,3,4)
CELLS=((F(0),F(1,32)),(F(7,32),F(1,4)),(F(15,32),F(1,2)),(F(31,32),F(1)))

def ivrat(q):
 q=F(q); return IV.mpf(q.numerator)/q.denominator

def bounds(v):
 # mpmath's documented internal binary tuples yield exact rational endpoints.
 from mpmath.libmp import to_rational
 return tuple(F(*to_rational(z)) for z in v._mpi_)

def conv(a,b,k): return sum((a[i]*b[k-i] for i in range(k+1)),IV.mpf(0))

def tanhjet(z):
 y=[1-2/(IV.exp(2*z[0])+1)]
 # y'=z'*(1-y*y): coefficient n only needs already computed y[:n+1].
 for n in range(len(z)-1):
  result=IV.mpf(0)
  for k in range(n+1):
   m=n-k; c=(IV.mpf(1) if m==0 else IV.mpf(0))-conv(y,y,m)
   result+=(k+1)*z[k+1]*c
  y.append(result/(n+1))
 return y

def interval_reference(net,a,b):
 x=IV.mpf([ivrat(a).a,ivrat(b).b]); order=6
 state=[[x,IV.mpf(1)]+[IV.mpf(0)]*(order-1)]
 for layer in net['layers']:
  result=[]
  for weights,bias in zip(layer['weight'],layer['bias']):
   v=[ivrat(F.from_float(float(bias)))]+[IV.mpf(0)]*order
   for w,jet in zip(weights,state):
    w=ivrat(F.from_float(float(w))); v=[v[j]+w*jet[j] for j in range(order+1)]
   result.append(tanhjet(v) if layer['activation']=='tanh' else v)
  state=result
 n=state[0]; q=[x*(1-x),1-2*x,IV.mpf(-1)]+[IV.mpf(0)]*4
 ans={}
 for j in ORDERS:
  val=-factorial(j+2)*conv(q,n,j+2)
  for i,c in enumerate(net['source_coefficients']):
   if i>=j:val-=ivrat(F.from_float(float(c)))*(factorial(i)//factorial(i-j))*x**(i-j)
  ans[j]=bounds(val)
 return ans

def atomic(path,obj):
 data=(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n').encode()
 tmp=path.with_suffix('.tmp')
 with tmp.open('xb') as h:h.write(data);h.flush();os.fsync(h.fileno())
 os.replace(tmp,path)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--fixtures',type=Path,nargs='+',required=True);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
 args.out.mkdir(parents=True,exist_ok=False)
 source=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0]
 needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
 if source.count(needle)!=1:raise ValueError('prototype structure changed')
 ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(source.replace(needle,'OUT=None'),'polydefs','exec'),ns)
 spec=importlib.util.spec_from_file_location('point_reference',ROOT/'scripts/107_polynomial_derivative_crosscheck_20261101.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
 points=sorted({x for a,b in CELLS for x in (a,(a+b)/2,b)});rows=[];fixtures=[];intervals=[]
 for path in args.fixtures:
  raw=path.read_bytes();sha=hashlib.sha256(raw).hexdigest();nets=json.loads(raw)['networks']; fixtures.append({'sha256':sha,'networks':len(nets)})
  if len({n['label'] for n in nets})!=len(nets):raise ValueError('duplicate labels')
  for net in nets:
   p=ns['PolynomialForcing'](path,net['label'],sha,net['source_coefficients']); fns={0:p.residual,2:p.r_second,3:p.r_third,4:p.r_fourth}
   for a,b in CELLS:
    cell=point(a).lower().union(point(b).upper()); indep=interval_reference(net,a,b)
    for j in ORDERS:
     ball=fns[j](cell);lo,hi=exact_fraction(ball.lower()),exact_fraction(ball.upper()); il,iu=indep[j]
     overlap=lo<=iu and il<=hi
     intervals.append({'fixture_sha256':sha,'network':net['label'],'cell':[str(a),str(b)],'order':j,'arb_lower':str(lo),'arb_upper':str(hi),'mpiv_lower':str(il),'mpiv_upper':str(iu),'overlap':overlap})
     for x in (a,(a+b)/2,b):
      low=mod.mp_value(net,x,j,60);high=mod.mp_value(net,x,j,90)
      with mp.workdps(90):
       stable=abs(low-high)<=mp.mpf('1e-50')*max(1,abs(high));value=mp.mpf(high)
       included=mp.mpf(lo.numerator)/lo.denominator<=value<=mp.mpf(hi.numerator)/hi.denominator
       included_iv=mp.mpf(il.numerator)/il.denominator<=value<=mp.mpf(iu.numerator)/iu.denominator
      rows.append({'fixture_sha256':sha,'network':net['label'],'x':str(x),'order':j,'cell':[str(a),str(b)],'mp_60_vs_90_stable':stable,'mp90_contained':included,'mpiv90_contained':included_iv})
   print(net['label'],'done',flush=True)
 atomic(args.out/'records.json',rows);atomic(args.out/'intervals.json',intervals)
 expected=sum(a['networks'] for a in fixtures)*len(points)*len(ORDERS); expected_cells=sum(a['networks'] for a in fixtures)*len(CELLS)*len(ORDERS)
 if len(rows)!=expected or len(intervals)!=expected_cells:raise ValueError('unexpected row count')
 checks={'fixtures':fixtures,'points':[str(x) for x in points],'orders':ORDERS,'expected':expected,'records':len(rows),'stable':sum(x['mp_60_vs_90_stable'] for x in rows),'contained':sum(x['mp90_contained'] for x in rows),'mpiv_contained':sum(x['mpiv90_contained'] for x in rows),'interval_expected':expected_cells,'interval_records':len(intervals),'interval_overlap':sum(x['overlap'] for x in intervals),'records_sha256':hashlib.sha256((args.out/'records.json').read_bytes()).hexdigest(),'intervals_sha256':hashlib.sha256((args.out/'intervals.json').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'failures':[r for r in rows if not (r['mp_60_vs_90_stable'] and r['mp90_contained'] and r['mpiv90_contained'])]+[r for r in intervals if not r['overlap']],'scope':'full-cell interval comparison and sampled numerical diagnostics; overlap or sampling does not prove either backend'}
 atomic(args.out/'summary.json',checks);print(json.dumps(checks),flush=True);return int(bool(checks['failures']))
if __name__=='__main__':sys.exit(main())
