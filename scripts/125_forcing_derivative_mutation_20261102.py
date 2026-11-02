"""Post-hoc sensitivity: omit f^(j) from stored residual derivative bounds.

For nonnegative monomial source coefficients on [0,1], analytically bound f^(j)
on each cell; add this to the recorded correct residual bound. This is a
conservative enclosure of the simulated omitted-forcing implementation, not
an actual execution of the verifier with mutated code.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
from collections import Counter
import argparse,hashlib,json,os,sys
CELLS=((F(0),F(1,32)),(F(7,32),F(1,4)),(F(15,32),F(1,2)),(F(31,32),F(1)))
SOURCES={'quadratic_12x2':(12,2),'cubic_20x3':(20,3),'quintic_42x5':(42,5)}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 data=(json.dumps(obj,indent=2)+'\n').encode();tmp=p.with_suffix('.tmp')
 with tmp.open('xb') as f:f.write(data);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 s=json.loads((a.input/'summary.json').read_text());rp=a.input/'records.json'
 if sha(rp)!=s['records_sha256']:raise ValueError('input hash mismatch')
 rows=[]
 for r in json.loads(rp.read_text()):
  left,right=CELLS[r['cell_index']];c,p=SOURCES[r['source']];j=r['order'];mult=F(c*factorial(p),factorial(p-j)) if p>=j else F(0);fl,fu=(mult*left**(p-j),mult*right**(p-j)) if p>=j else (F(0),F(0));lo,hi=F(r['arb_lower'])+fl,F(r['arb_upper'])+fu;point=F(r['mp70']);detected=not lo<=point<=hi
  rows.append({'fixture_sha256':r['fixture_sha256'],'network':r['network'],'source':r['source'],'cell_index':r['cell_index'],'order':j,'mutated_lower':str(lo),'mutated_upper':str(hi),'point':str(point),'forcing_derivative_nonzero':mult!=0,'detected':detected})
 atom(a.out/'records.json',rows);group={}
 for source in SOURCES:
  for order in (0,2,3,4):
   subset=[r for r in rows if r['source']==source and r['order']==order];group[source+':'+str(order)]={'detected':sum(r['detected'] for r in subset),'total':len(subset),'nonzero':sum(r['forcing_derivative_nonzero'] for r in subset)}
 summary={'expected':s['records'],'records':len(rows),'detected':sum(r['detected'] for r in rows),'nonzero':sum(r['forcing_derivative_nonzero'] for r in rows),'by_source_order':group,'input_sha256':s['records_sha256'],'records_sha256':sha(a.out/'records.json'),'script_sha256':sha(Path(__file__)),'scope':'post-hoc conservative simulated omission; numerical point references only; a missed detection does not establish correctness'}
 if len(rows)!=s['records']:raise ValueError('count')
 atom(a.out/'summary.json',summary);print(json.dumps(summary));return 0
if __name__=='__main__':sys.exit(main())
