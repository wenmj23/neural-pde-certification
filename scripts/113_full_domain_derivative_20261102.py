"""Pre-specified full-domain crosscheck of interval residual derivatives.

32 equal cells, 24 existing frozen networks, derivative orders 0,2,3,4.
Two interval propagations share the PDE and inputs; midpoint checks are numerical.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,importlib.util,json,os,sys
import mpmath as mp
from arb_certification import exact_fraction
from arb_reconstruction import point
ROOT=Path(__file__).resolve().parents[1]
def module(name,file):
 s=importlib.util.spec_from_file_location(name,ROOT/'scripts'/file);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
iv=module('ivcheck','111_interval_derivative_crosscheck_20261102.py')
pointref=module('pointcheck','107_polynomial_derivative_crosscheck_20261101.py')
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def atom(path,obj):
 b=(json.dumps(obj,indent=2,allow_nan=False)+'\n').encode();tmp=path.with_suffix('.tmp')
 with tmp.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(tmp,path)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--fixtures',nargs='+',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 source=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
 if source.count(needle)!=1:raise ValueError('prototype changed')
 ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(source.replace(needle,'OUT=None'),'polydefs','exec'),ns)
 rows=[];fixtures=[]
 for path in a.fixtures:
  digest=sha(path);nets=json.loads(path.read_bytes())['networks'];fixtures.append({'sha256':digest,'networks':len(nets)})
  if len({n['label'] for n in nets})!=len(nets):raise ValueError('duplicate labels')
  for net in nets:
   p=ns['PolynomialForcing'](path,net['label'],digest,net['source_coefficients']); fns={0:p.residual,2:p.r_second,3:p.r_third,4:p.r_fourth}
   for k in range(32):
    left,right=F(k,32),F(k+1,32);mid=(left+right)/2;cell=point(left).lower().union(point(right).upper());independent=iv.interval_reference(net,left,right)
    for order in (0,2,3,4):
     v=fns[order](cell);lo,hi=exact_fraction(v.lower()),exact_fraction(v.upper());il,iu=independent[order]
     low=pointref.mp_value(net,mid,order,50);high=pointref.mp_value(net,mid,order,70)
     with mp.workdps(70):
      stable=abs(low-high)<=mp.mpf('1e-40')*max(1,abs(high));x=mp.mpf(high)
      inarb=mp.mpf(lo.numerator)/lo.denominator<=x<=mp.mpf(hi.numerator)/hi.denominator
      iniv=mp.mpf(il.numerator)/il.denominator<=x<=mp.mpf(iu.numerator)/iu.denominator
     rows.append({'fixture_sha256':digest,'network':net['label'],'cell_index':k,'order':order,'arb_lower':str(lo),'arb_upper':str(hi),'mpiv_lower':str(il),'mpiv_upper':str(iu),'interval_overlap':lo<=iu and il<=hi,'midpoint_stable':stable,'midpoint_in_arb':inarb,'midpoint_in_mpiv':iniv})
   print(net['label'],'32 cells completed',flush=True)
 atom(a.out/'records.json',rows);expected=sum(f['networks'] for f in fixtures)*32*4
 if len(rows)!=expected:raise ValueError('unexpected count')
 s={'fixtures':fixtures,'cells':32,'orders':[0,2,3,4],'expected':expected,'records':len(rows),'interval_overlap':sum(r['interval_overlap'] for r in rows),'midpoint_stable':sum(r['midpoint_stable'] for r in rows),'midpoint_in_arb':sum(r['midpoint_in_arb'] for r in rows),'midpoint_in_mpiv':sum(r['midpoint_in_mpiv'] for r in rows),'records_sha256':sha(a.out/'records.json'),'script_sha256':sha(Path(__file__)),'independent_interval_script_sha256':sha(ROOT/'scripts/111_interval_derivative_crosscheck_20261102.py'),'point_reference_script_sha256':sha(ROOT/'scripts/107_polynomial_derivative_crosscheck_20261101.py'),'failures':[{'network':r['network'],'cell_index':r['cell_index'],'order':r['order']} for r in rows if not all(r[x] for x in ('interval_overlap','midpoint_stable','midpoint_in_arb','midpoint_in_mpiv'))],'scope':'finite case diagnostics, not a proof of interval backends'}
 atom(a.out/'summary.json',s);print(json.dumps(s),flush=True);return int(bool(s['failures']))
if __name__=='__main__':sys.exit(main())
