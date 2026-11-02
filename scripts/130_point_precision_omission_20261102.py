"""Post-hoc precision and exact-midpoint check for the remaining edge cases."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,os,sys
from flint import arb
from arb_certification import exact_fraction
from arb_reconstruction import point
ROOT=Path(__file__).resolve().parents[1]
SOURCES={'quadratic_12x2':[0,0,12],'cubic_20x3':[0,0,0,20],'quintic_42x5':[0,0,0,0,0,42]};CELLS=((F(0),F(1,32)),(F(7,32),F(1,4)),(F(15,32),F(1,2)),(F(31,32),F(1)))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2)+'\n').encode();t=p.with_suffix('.tmp')
 with t.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(t,p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True,type=Path);ap.add_argument('--fixtures',nargs='+',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 s=json.loads((a.input/'summary.json').read_text());p=a.input/'records.json'
 if sha(p)!=s['records_sha256']:raise ValueError('parent hash mismatch')
 selected=[x for x in json.loads(p.read_text()) if x['shrink']==4096 and not x['mutation_detected']];fixtures={sha(p):p for p in a.fixtures}
 if any(x['fixture_sha256'] not in fixtures for x in selected):raise ValueError('fixture mismatch')
 src=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
 if src.count(needle)!=1:raise ValueError('prototype changed')
 ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(src.replace(needle,'OUT=None'),'defs','exec'),ns)
 class Omit(ns['PolynomialForcing']):
  def source_derivative(self,x,j):return arb(0)
 rows=[];orders={0:'residual',2:'r_second',3:'r_third',4:'r_fourth'}
 for r in selected:
  h,label,source=r['fixture_sha256'],r['network'],r['source'];mid=sum(CELLS[r['cell_index']],F(0))/2;reference=F(r['midpoint_reference'])
  for bits in (103,160):
   valid=ns['PolynomialForcing'](fixtures[h],label,h,SOURCES[source],bits);broken=Omit(fixtures[h],label,h,SOURCES[source],bits)
   cell=point(mid);v=getattr(valid,orders[r['order']])(cell);b=getattr(broken,orders[r['order']])(cell);vl,vu=exact_fraction(v.lower()),exact_fraction(v.upper());bl,bu=exact_fraction(b.lower()),exact_fraction(b.upper())
   rows.append({'fixture_sha256':h,'network':label,'source':source,'cell_index':r['cell_index'],'order':r['order'],'bits':bits,'correct_lower':str(vl),'correct_upper':str(vu),'mutated_lower':str(bl),'mutated_upper':str(bu),'midpoint_reference':str(reference),'correct_contains':vl<=reference<=vu,'mutation_detected':not bl<=reference<=bu})
 atom(a.out/'records.json',rows);summary={'input_sha256':s['records_sha256'],'selected_count':len(selected),'bits':[103,160],'expected':len(selected)*2,'records':len(rows),'by_bits':{str(i):{'total':sum(r['bits']==i for r in rows),'correct_contains':sum(r['correct_contains'] for r in rows if r['bits']==i),'mutation_detected':sum(r['mutation_detected'] for r in rows if r['bits']==i)} for i in (103,160)},'records_sha256':sha(a.out/'records.json'),'script_sha256':sha(Path(__file__)),'failures':[r for r in rows if not r['correct_contains']],'scope':'post-hoc exact point, not interval coverage; reference is high-precision numeric approximation'}
 if len(rows)!=summary['expected']:raise ValueError('count mismatch')
 atom(a.out/'summary.json',summary);print(json.dumps(summary));return int(bool(summary['failures']))
if __name__=='__main__':sys.exit(main())
