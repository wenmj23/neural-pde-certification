"""Execute an intentionally faulty source-omission subclass on nested cells.

The production code is untouched. Reference midpoint values are from 122;
interval bounds are recomputed on intervals centered at those same points.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,os,sys
from flint import arb
from arb_certification import exact_fraction
from arb_reconstruction import point
ROOT=Path(__file__).resolve().parents[1]
CELLS=((F(0),F(1,32)),(F(7,32),F(1,4)),(F(15,32),F(1,2)),(F(31,32),F(1)))
POWER={'quadratic_12x2':2,'cubic_20x3':3,'quintic_42x5':5}
SOURCES={'quadratic_12x2':[0,0,12],'cubic_20x3':[0,0,0,20],'quintic_42x5':[0,0,0,0,0,42]}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2)+'\n').encode();tmp=p.with_suffix('.tmp')
 with tmp.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True,type=Path);ap.add_argument('--fixtures',required=True,nargs='+',type=Path);ap.add_argument('--out',required=True,type=Path);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 s=json.loads((a.input/'summary.json').read_text());rpath=a.input/'records.json'
 if sha(rpath)!=s['records_sha256']:raise ValueError('input hash mismatch')
 fixt={sha(p):p for p in a.fixtures}
 if set(fixt)!={x['sha256'] for x in s['fixtures']}:raise ValueError('fixture hash mismatch')
 src=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
 if src.count(needle)!=1:raise ValueError('prototype layout changed')
 ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(src.replace(needle,'OUT=None'),'defs','exec'),ns)
 class OmitSource(ns['PolynomialForcing']):
  def source_derivative(self,x,j):return arb(0)
 rows=[];cache={};orders={0:'residual',2:'r_second',3:'r_third',4:'r_fourth'}
 for rec in json.loads(rpath.read_text()):
  order=rec['order'];source=rec['source']
  if order>POWER[source]:continue
  identity=(rec['fixture_sha256'],rec['network'],source)
  if identity not in cache:
   h,label,source=identity;cache[identity]=(ns['PolynomialForcing'](fixt[h],label,h,SOURCES[source]),OmitSource(fixt[h],label,h,SOURCES[source]))
  valid,broken=cache[identity];left,right=CELLS[rec['cell_index']];mid=(left+right)/2;ref=F(rec['mp70'])
  for shrink in (1,2,4,8,16):
   half=(right-left)/2/shrink;lo,hi=mid-half,mid+half;cell=point(lo).lower().union(point(hi).upper())
   v=getattr(valid,orders[order])(cell);b=getattr(broken,orders[order])(cell);vl,vu=exact_fraction(v.lower()),exact_fraction(v.upper());bl,bu=exact_fraction(b.lower()),exact_fraction(b.upper())
   rows.append({'fixture_sha256':identity[0],'network':identity[1],'source':source,'cell_index':rec['cell_index'],'order':order,'shrink':shrink,'correct_lower':str(vl),'correct_upper':str(vu),'mutated_lower':str(bl),'mutated_upper':str(bu),'midpoint_reference':str(ref),'correct_contains':vl<=ref<=vu,'mutation_detected':not bl<=ref<=bu})
  if len(rows)%720==0:print('completed',len(rows),'of 4320 checks',flush=True)
 atom(a.out/'records.json',rows);summary={'input_sha256':s['records_sha256'],'fixtures':s['fixtures'],'shrink':[1,2,4,8,16],'expected':864*5,'records':len(rows),'by_shrink':{str(i):{'total':sum(r['shrink']==i for r in rows),'correct_contains':sum(r['correct_contains'] for r in rows if r['shrink']==i),'mutation_detected':sum(r['mutation_detected'] for r in rows if r['shrink']==i)} for i in (1,2,4,8,16)},'records_sha256':sha(a.out/'records.json'),'script_sha256':sha(Path(__file__)),'failures':[{'network':r['network'],'source':r['source'],'cell_index':r['cell_index'],'order':r['order'],'shrink':r['shrink']} for r in rows if not r['correct_contains']],'scope':'intentional fault injection, nested cells centered on prior numerical references; diagnostic not formal proof'}
 if len(rows)!=summary['expected']:raise ValueError('unexpected count')
 atom(a.out/'summary.json',summary);print(json.dumps(summary),flush=True);return int(bool(summary['failures']))
if __name__=='__main__':sys.exit(main())
