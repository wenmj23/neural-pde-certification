"""Test direct point-evaluation check on all 864 nonzero-source cases.
Faulty subclass is isolated; saved 70-digit values are numerical references.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,os,sys
from flint import arb
from arb_certification import exact_fraction
from arb_reconstruction import point
ROOT=Path(__file__).resolve().parents[1];CELLS=((F(0),F(1,32)),(F(7,32),F(1,4)),(F(15,32),F(1,2)),(F(31,32),F(1)));SOURCES={'quadratic_12x2':[0,0,12],'cubic_20x3':[0,0,0,20],'quintic_42x5':[0,0,0,0,0,42]};POWERS={'quadratic_12x2':2,'cubic_20x3':3,'quintic_42x5':5}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2)+'\n').encode();t=p.with_suffix('.tmp')
 with t.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(t,p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True,type=Path);ap.add_argument('--fixtures',required=True,nargs='+',type=Path);ap.add_argument('--out',required=True,type=Path);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 s=json.loads((a.input/'summary.json').read_text());ip=a.input/'records.json'
 if sha(ip)!=s['records_sha256']:raise ValueError('input hash')
 fixtures={sha(p):p for p in a.fixtures}
 if set(fixtures)!={f['sha256'] for f in s['fixtures']}:raise ValueError('fixture mismatch')
 src=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
 if src.count(needle)!=1:raise ValueError('prototype changed')
 ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(src.replace(needle,'OUT=None'),'defs','exec'),ns)
 class Omit(ns['PolynomialForcing']):
  def source_derivative(self,x,j):return arb(0)
 orders={0:'residual',2:'r_second',3:'r_third',4:'r_fourth'};cache={};rows=[]
 for r in json.loads(ip.read_text()):
  h,label,source=r['fixture_sha256'],r['network'],r['source'];j=r['order']
  if j>POWERS[source]:continue
  k=(h,label,source)
  if k not in cache:cache[k]=(ns['PolynomialForcing'](fixtures[h],label,h,SOURCES[source]),Omit(fixtures[h],label,h,SOURCES[source]))
  valid,broken=cache[k];mid=sum(CELLS[r['cell_index']],F(0))/2;ref=F(r['mp70']);v=getattr(valid,orders[j])(point(mid));b=getattr(broken,orders[j])(point(mid));vl,vu=exact_fraction(v.lower()),exact_fraction(v.upper());bl,bu=exact_fraction(b.lower()),exact_fraction(b.upper());rows.append({'fixture_sha256':h,'network':label,'source':source,'cell_index':r['cell_index'],'order':j,'correct_lower':str(vl),'correct_upper':str(vu),'mutated_lower':str(bl),'mutated_upper':str(bu),'midpoint_reference':str(ref),'correct_contains':vl<=ref<=vu,'mutation_detected':not bl<=ref<=bu})
 atom(a.out/'records.json',rows);s={'input_sha256':s['records_sha256'],'expected':864,'records':len(rows),'correct_contains':sum(r['correct_contains'] for r in rows),'mutation_detected':sum(r['mutation_detected'] for r in rows),'records_sha256':sha(a.out/'records.json'),'script_sha256':sha(Path(__file__)),'failures':[r for r in rows if not r['correct_contains']],'scope':'post-hoc exact-point, simulated source omission, numerical point reference; no assertion for unsampled x'}
 if len(rows)!=s['expected']:raise ValueError('count')
 atom(a.out/'summary.json',s);print(json.dumps(s));return int(bool(s['failures']))
if __name__=='__main__':sys.exit(main())
