"""Cross-check residual derivatives for transferred polynomial forcings.
4 fixed cells, 4 derivative orders; save exact numerical midpoint references.
"""
from pathlib import Path
from fractions import Fraction as F
from copy import deepcopy
import argparse,hashlib,importlib.util,json,os,sys
from mpmath.libmp import to_rational
from arb_certification import exact_fraction
from arb_reconstruction import point
ROOT=Path(__file__).resolve().parents[1]
SOURCES={'quadratic_12x2':[0,0,12],'cubic_20x3':[0,0,0,20],'quintic_42x5':[0,0,0,0,0,42]};CELLS=((F(0),F(1,32)),(F(7,32),F(1,4)),(F(15,32),F(1,2)),(F(31,32),F(1)))
def mod(name,file):
 s=importlib.util.spec_from_file_location(name,ROOT/'scripts'/file);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
iv=mod('ivcheck','111_interval_derivative_crosscheck_20261102.py');mpref=mod('mpref','107_polynomial_derivative_crosscheck_20261101.py')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2)+'\n').encode();tmp=p.with_suffix('.tmp')
 with tmp.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--fixtures',type=Path,nargs='+',required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 src=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
 if src.count(needle)!=1:raise ValueError('prototype changed')
 ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(src.replace(needle,'OUT=None'),'defs','exec'),ns)
 fixtures=[];rows=[]
 for path in a.fixtures:
  digest=sha(path);nets=json.loads(path.read_bytes())['networks'];fixtures.append({'sha256':digest,'networks':len(nets)})
  for net in nets:
   for source,coeff in SOURCES.items():
    p=ns['PolynomialForcing'](path,net['label'],digest,coeff);fns={0:p.residual,2:p.r_second,3:p.r_third,4:p.r_fourth};raw=deepcopy(net);raw['source_coefficients']=coeff
    for i,(left,right) in enumerate(CELLS):
     cell=point(left).lower().union(point(right).upper());indep=iv.interval_reference(raw,left,right);mid=(left+right)/2
     for order in (0,2,3,4):
      ball=fns[order](cell);lo,hi=exact_fraction(ball.lower()),exact_fraction(ball.upper());il,iu=indep[order]
      low=F(*to_rational(mpref.mp_value(raw,mid,order,50)._mpf_));high=F(*to_rational(mpref.mp_value(raw,mid,order,70)._mpf_));stable=abs(low-high)<=F(1,10**40)*max(F(1),abs(high))
      rows.append({'fixture_sha256':digest,'network':net['label'],'source':source,'cell_index':i,'order':order,'arb_lower':str(lo),'arb_upper':str(hi),'mpiv_lower':str(il),'mpiv_upper':str(iu),'mp50':str(low),'mp70':str(high),'overlap':lo<=iu and il<=hi,'stable':stable,'point_in_arb':lo<=high<=hi,'point_in_mpiv':il<=high<=iu})
   print(net['label'],'done',flush=True)
 atom(a.out/'records.json',rows);n=sum(x['networks'] for x in fixtures);expected=n*len(SOURCES)*len(CELLS)*4
 s={'fixtures':fixtures,'sources':SOURCES,'cells':[[str(a),str(b)] for a,b in CELLS],'orders':[0,2,3,4],'expected':expected,'records':len(rows),'overlap':sum(r['overlap'] for r in rows),'stable':sum(r['stable'] for r in rows),'point_in_arb':sum(r['point_in_arb'] for r in rows),'point_in_mpiv':sum(r['point_in_mpiv'] for r in rows),'records_sha256':sha(a.out/'records.json'),'script_sha256':sha(Path(__file__)),'failures':[{'network':r['network'],'source':r['source'],'cell':r['cell_index'],'order':r['order']} for r in rows if not all(r[k] for k in ('overlap','stable','point_in_arb','point_in_mpiv'))],'scope':'finite interval and point cross-check; mpmath point derivatives are numerical approximations, not independently certified'}
 if len(rows)!=expected:raise ValueError('incorrect count')
 atom(a.out/'summary.json',s);print(json.dumps(s),flush=True);return int(bool(s['failures']))
if __name__=='__main__':sys.exit(main())
