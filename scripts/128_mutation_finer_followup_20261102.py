"""Post-hoc investigate faults not detected at 1/16 original width."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,os,sys
from flint import arb
from arb_certification import exact_fraction
from arb_reconstruction import point
ROOT=Path(__file__).resolve().parents[1]; CELLS=((F(0),F(1,32)),(F(7,32),F(1,4)),(F(15,32),F(1,2)),(F(31,32),F(1))); SOURCES={'quadratic_12x2':[0,0,12],'cubic_20x3':[0,0,0,20],'quintic_42x5':[0,0,0,0,0,42]}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2)+'\n').encode();tmp=p.with_suffix('.tmp')
 with tmp.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True,type=Path);ap.add_argument('--fixtures',required=True,nargs='+',type=Path);ap.add_argument('--out',required=True,type=Path);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 s=json.loads((a.input/'summary.json').read_text());p=a.input/'records.json'
 if sha(p)!=s['records_sha256']:raise ValueError('parent hash mismatch')
 selected=[x for x in json.loads(p.read_text()) if x['shrink']==16 and not x['mutation_detected']];fixt={sha(p):p for p in a.fixtures}
 if set(fixt)!={f['sha256'] for f in s['fixtures']}:raise ValueError('fixtures differ')
 src=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
 if src.count(needle)!=1:raise ValueError('prototype changed')
 ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(src.replace(needle,'OUT=None'),'defs','exec'),ns)
 class Omit(ns['PolynomialForcing']):
  def source_derivative(self,x,j):return arb(0)
 cache={};orders={0:'residual',2:'r_second',3:'r_third',4:'r_fourth'};rows=[]
 for r in selected:
  key=(r['fixture_sha256'],r['network'],r['source']);h,label,source=key
  if key not in cache:cache[key]=(ns['PolynomialForcing'](fixt[h],label,h,SOURCES[source]),Omit(fixt[h],label,h,SOURCES[source]))
  valid,broken=cache[key];left,right=CELLS[r['cell_index']];mid=(left+right)/2;ref=F(r['midpoint_reference'])
  for shrink in (32,64,128,256,1024,4096):
   half=(right-left)/2/shrink;cell=point(mid-half).lower().union(point(mid+half).upper());v=getattr(valid,orders[r['order']])(cell);b=getattr(broken,orders[r['order']])(cell);vl,vu=exact_fraction(v.lower()),exact_fraction(v.upper());bl,bu=exact_fraction(b.lower()),exact_fraction(b.upper());rows.append({'fixture_sha256':h,'network':label,'source':source,'cell_index':r['cell_index'],'order':r['order'],'shrink':shrink,'correct_lower':str(vl),'correct_upper':str(vu),'mutated_lower':str(bl),'mutated_upper':str(bu),'midpoint_reference':str(ref),'correct_contains':vl<=ref<=vu,'mutation_detected':not bl<=ref<=bu})
 atom(a.out/'records.json',rows);summary={'input_sha256':s['records_sha256'],'selected_count':len(selected),'shrink':[32,64,128,256,1024,4096],'expected':len(selected)*6,'records':len(rows),'by_shrink':{str(i):{'total':sum(r['shrink']==i for r in rows),'correct_contains':sum(r['correct_contains'] for r in rows if r['shrink']==i),'mutation_detected':sum(r['mutation_detected'] for r in rows if r['shrink']==i)} for i in (32,64,128,256,1024,4096)},'records_sha256':sha(a.out/'records.json'),'script_sha256':sha(Path(__file__)),'failures':[r for r in rows if not r['correct_contains']],'scope':'post-hoc selected failures at factor 16; numerical point diagnostic, not general detection rate'}
 if len(rows)!=summary['expected']:raise ValueError('count mismatch')
 atom(a.out/'summary.json',summary);print(json.dumps(summary),flush=True);return int(bool(summary['failures']))
if __name__=='__main__':sys.exit(main())
