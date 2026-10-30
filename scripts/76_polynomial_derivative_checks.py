"""High precision point diagnostics of new residual and source derivatives."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import json,hashlib
import mpmath as mp
from arb_reconstruction import point
from arb_certification import exact_fraction
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261030/polynomial_derivative_checks';OUT.mkdir(parents=True,exist_ok=False)
s=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];s=s.replace("OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)","OUT=None");ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(s,'poly_definitions','exec'),ns)
s=(ROOT/'scripts/62_derivative_crosscheck.py').read_text().split('for name in CASES:')[0];s=s.replace("OUT=ROOT.parent/'run_20261030/derivative_crosscheck';OUT.mkdir(parents=True,exist_ok=False)","OUT=None");nn={'__file__':str(ROOT/'scripts/62_derivative_crosscheck.py')};exec(compile(s,'mp_network_definitions','exec'),nn);make_fn=nn['make_fn'];mp.mp.dps=90
fixture=ROOT/'configs/frozen_networks_holdout_2026-10-30.json';data=fixture.read_bytes();sha=hashlib.sha256(data).hexdigest();nets=json.loads(data)['networks'];rows=[]
def mpq(q):return mp.mpf(q.numerator)/q.denominator
def fder(coeff,x,j):return sum((mp.mpf(coeff[i]*(factorial(i)//factorial(i-j)))*x**(i-j) for i in range(j,len(coeff))),mp.mpf(0))
for net in nets:
 fn=make_fn(net);samples={}
 for center in (F(1,7),F(1,2),F(6,7)):
  a,b=center-F(1,1000),center+F(1,1000)
  for q in (a,center,b):samples[q]=mp.taylor(fn,mpq(q),6)
 for name,coeff in ns['SOURCES'].items():
  prob=ns['PolynomialForcing'](fixture,net['label'],sha,coeff)
  for center in (F(1,7),F(1,2),F(6,7)):
   a,b=center-F(1,1000),center+F(1,1000);x=point(a).lower().union(point(b).upper());ivs={0:prob.residual(x),2:prob.r_second(x),3:prob.r_third(x),4:prob.r_fourth(x)}
   for q in (a,center,b):
    for j in (0,2,3,4):
     val=-factorial(j+2)*samples[q][j+2]-fder(coeff,mpq(q),j);lo=exact_fraction(ivs[j].lower());hi=exact_fraction(ivs[j].upper());rows.append({'network':net['label'],'source':name,'interval_center':str(center),'point':str(q),'order':j,'lower':str(lo),'upper':str(hi),'value':mp.nstr(val,80),'passed':mpq(lo)<=val<=mpq(hi)})
  (OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n');print(net['label'],name,'done',flush=True)
summary={'networks':len(nets),'sources':len(ns['SOURCES']),'checks':len(rows),'passed':sum(x['passed'] for x in rows),'failures':[x for x in rows if not x['passed']],'caveat':'90-digit mpmath point diagnostics, not an interval proof','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k!='failures'}),'failures',len(summary['failures']),flush=True)
