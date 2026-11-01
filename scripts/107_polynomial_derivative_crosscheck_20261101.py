"""Cross-check polynomial residual derivatives against an independent mp network.

Numerical mp.diff at two precisions is a diagnostic, not a rigorous enclosure.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import argparse,hashlib,json,os,sys
import mpmath as mp
from arb_reconstruction import point
from arb_certification import exact_fraction
ROOT=Path(__file__).resolve().parents[1]
def atom(p,obj):
 b=(json.dumps(obj,indent=2)+'\n').encode();tmp=p.with_suffix(p.suffix+'.tmp')
 with tmp.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def mp_value(net,rat,order,precision):
 with mp.workdps(precision):
  layers=[([[mp.mpf(x) for x in row] for row in layer['weight']], [mp.mpf(x) for x in layer['bias']],layer['activation']) for layer in net['layers']]
  def solution(t):
   a=[t]
   for weights,bias,activation in layers:
    a=[sum((w*z for w,z in zip(row,a)),b) for row,b in zip(weights,bias)]
    if activation=='tanh':a=[mp.tanh(z) for z in a]
    elif activation!='identity':raise ValueError('Unsupported activation')
   if len(a)!=1:raise ValueError('non-scalar network')
   return t*(1-t)*a[0]
  x=mp.mpf(rat.numerator)/rat.denominator
  forcing=sum((mp.mpf(c)*factorial(i)/factorial(i-order)*x**(i-order) for i,c in enumerate(net['source_coefficients']) if i>=order),mp.mpf(0))
  return +(-mp.diff(solution,x,order+2)-forcing)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--fixtures',nargs='+',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
 src=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];target="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
 if src.count(target)!=1:raise ValueError('prototype changed')
 ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(src.replace(target,'OUT=None'),'polydefs','exec'),ns)
 grid=[F(1,7),F(1,2),F(6,7)];orders=(0,2,3,4);rows=[];fixtures=[]
 for path in args.fixtures:
  path=path.resolve();raw=path.read_bytes();sha=hashlib.sha256(raw).hexdigest();nets=json.loads(raw)['networks'];fixtures.append({'sha256':sha,'networks':len(nets)})
  for net in nets:
   p=ns['PolynomialForcing'](path,net['label'],sha,net['source_coefficients']);fns={0:p.residual,2:p.r_second,3:p.r_third,4:p.r_fourth}
   for x in grid:
    for order in orders:
     ball=fns[order](point(x));lo,hi=exact_fraction(ball.lower()),exact_fraction(ball.upper());low=mp_value(net,x,order,60);high=mp_value(net,x,order,90)
     with mp.workdps(90):
      left=mp.mpf(lo.numerator)/lo.denominator;right=mp.mpf(hi.numerator)/hi.denominator;stable=abs(low-high)<=mp.mpf('1e-50')*max(1,abs(high));contained=left<=high<=right
     rows.append({'fixture_sha256':sha,'network':net['label'],'x':str(x),'order':order,'arb_lower':str(lo),'arb_upper':str(hi),'mp_60_vs_90_stable':stable,'mp90_contained':contained})
   print(net['label'],'completed',flush=True)
 atom(out/'records.json',rows)
 saved=json.loads((out/'records.json').read_text()); expected=sum(a['networks'] for a in fixtures)*len(grid)*len(orders);keys={(r['fixture_sha256'],r['network'],r['x'],r['order']) for r in saved};assert len(saved)==len(keys)==expected
 checks={'fixtures':fixtures,'points':[str(x) for x in grid],'orders':orders,'expected':expected,'records':len(saved),'stable':sum(r['mp_60_vs_90_stable'] for r in saved),'contained':sum(r['mp90_contained'] for r in saved),'failures':[r for r in saved if not r['mp_60_vs_90_stable'] or not r['mp90_contained']],'records_sha256':hashlib.sha256((out/'records.json').read_bytes()).hexdigest(),'note':'mp.diff floating high precision diagnostic, not a rigorous reference or full interval proof'};atom(out/'summary.json',checks);print(json.dumps(checks),flush=True);return int(bool(checks['failures']))
if __name__=='__main__':sys.exit(main())
