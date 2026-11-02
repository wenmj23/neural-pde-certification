"""Pointwise sine residual derivative audit against independent mpmath.

Research diagnostic, not a proof on an interval. Fault injection is available
only through the imported run() API for the separate controlled-test scripts.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,os,sys
import mpmath as mp
from mpmath.libmp import to_rational
from arb_certification import exact_fraction,finite
from arb_cubic_reconstruction import SixthOrderDeep
from arb_reconstruction import point
POINTS=(F(1,7),F(2,7),F(5,7),F(1,11));ORDERS=(0,2,3,4)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2,allow_nan=False)+'\n').encode();tmp=p.with_suffix('.tmp')
 with tmp.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def source_zero(k,x,j):
 n,d=x.numerator,x.denominator
 return (k*n)%d==0 if j%2==0 else (2*k*n)%(2*d)==d
def mp_value(net,x,order,dps):
 with mp.workdps(dps):
  layers=[([[mp.mpf(v) for v in row] for row in layer['weight']],[mp.mpf(v) for v in layer['bias']],layer['activation']) for layer in net['layers']]
  def solution(t):
   state=[t]
   for weights,bias,activation in layers:
    state=[sum((w*z for w,z in zip(row,state)),b) for row,b in zip(weights,bias)]
    if activation=='tanh':state=[mp.tanh(z) for z in state]
    elif activation!='identity':raise ValueError('unsupported activation')
   if len(state)!=1:raise ValueError('non-scalar output')
   return t*(1-t)*state[0]
  xx=mp.mpf(x.numerator)/x.denominator;w=net['frequency']*mp.pi
  return +(-mp.diff(solution,xx,order+2)-w**(order+2)*mp.sin(w*xx+order*mp.pi/2))
class OmitForcing(SixthOrderDeep):
 """Deliberately faulty isolated subclass; all forcing terms omitted."""
 def residual(self,x):
  n=self.network_jet(x,2);return finite(-2*(self.q(x)*n[2]+(1-2*x)*n[1]-n[0]),'fault residual')
 def r_second(self,x):
  n=self.network_jet(x,4);return finite(-24*(self.q(x)*n[4]+(1-2*x)*n[3]-n[2]),'fault residual second')
 def r_third(self,x):
  n=self.network_jet(x,5);return finite(-120*(self.q(x)*n[5]+(1-2*x)*n[4]-n[3]),'fault residual third')
 def r_fourth(self,x):
  n=self.network_jet(x,6);return finite(-720*(self.q(x)*n[6]+(1-2*x)*n[5]-n[4]),'fault residual fourth')
def run(fixtures,out,omit_for_test=False):
 out.mkdir(parents=True,exist_ok=False);rows=[];inputs=[]
 for path in fixtures:
  h=sha(path);nets=json.loads(path.read_bytes())['networks'];labels={n['label'] for n in nets}
  if not labels or len(labels)!=len(nets):raise ValueError('empty or duplicate network labels')
  inputs.append({'sha256':h,'networks':len(nets)})
  for net in nets:
   k=net['frequency']
   if type(k) is not int or k<1:raise ValueError('invalid frequency')
   p=OmitForcing(path,net['label'],h) if omit_for_test else SixthOrderDeep(path,net['label'],h);fns={0:p.residual,2:p.r_second,3:p.r_third,4:p.r_fourth}
   for x in POINTS:
    for order in ORDERS:
     v=fns[order](point(x));lo,hi=exact_fraction(v.lower()),exact_fraction(v.upper());low=F(*to_rational(mp_value(net,x,order,60)._mpf_));high=F(*to_rational(mp_value(net,x,order,90)._mpf_));stable=abs(low-high)<=F(1,10**50)*max(F(1),abs(high));inside=lo<=high<=hi;zero=source_zero(k,x,order);expected_inside=not omit_for_test or zero
     rows.append({'fixture_sha256':h,'network':net['label'],'frequency':k,'x':str(x),'order':order,'arb_lower':str(lo),'arb_upper':str(hi),'mp60':str(low),'mp90':str(high),'source_zero_exact':zero,'stable':stable,'contained':inside,'expected_contained':expected_inside,'passed':stable and inside==expected_inside})
   print(net['label'],'completed',flush=True)
 atom(out/'records.json',rows);expected=sum(i['networks'] for i in inputs)*len(POINTS)*len(ORDERS);s={'inputs':inputs,'points':[str(p) for p in POINTS],'orders':ORDERS,'test_fault':'omit_forcing' if omit_for_test else 'none','expected':expected,'records':len(rows),'stable':sum(r['stable'] for r in rows),'contained':sum(r['contained'] for r in rows),'source_zero_exact':sum(r['source_zero_exact'] for r in rows),'expected_mismatch':sum(not r['expected_contained'] for r in rows),'passed':sum(r['passed'] for r in rows),'records_sha256':sha(out/'records.json'),'script_sha256':sha(Path(__file__)),'failures':[{'network':r['network'],'x':r['x'],'order':r['order']} for r in rows if not r['passed']],'scope':'finite points with numerical derivative reference; not validated continuous bound'}
 if len(rows)!=expected:raise ValueError('record count')
 atom(out/'summary.json',s);print(json.dumps(s),flush=True);return int(bool(s['failures']))
def main():
 p=argparse.ArgumentParser();p.add_argument('--fixtures',required=True,nargs='+',type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args();return run(a.fixtures,a.out)
if __name__=='__main__':sys.exit(main())
