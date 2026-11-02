"""Post-exit integrity and arithmetic verification for script 125."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import argparse,hashlib,json,sys
CELLS=[(F(i,32),F(i+1,32)) for i in (0,7,15,31)]
MONOMIALS={'quadratic_12x2':(12,2),'cubic_20x3':(20,3),'quintic_42x5':(42,5)}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True,type=Path);ap.add_argument('--mutation',required=True,type=Path);a=ap.parse_args();i,m=a.input,a.mutation;si=json.loads((i/'summary.json').read_text());sm=json.loads((m/'summary.json').read_text())
 if sha(i/'records.json')!=si['records_sha256'] or sha(m/'records.json')!=sm['records_sha256'] or sm['input_sha256']!=si['records_sha256']:raise ValueError('hash mismatch')
 inputs=json.loads((i/'records.json').read_text());outputs=json.loads((m/'records.json').read_text());key=lambda r:(r['fixture_sha256'],r['network'],r['source'],r['cell_index'],r['order']);old={key(r):r for r in inputs};out={key(r):r for r in outputs}
 if len(old)!=len(inputs) or len(out)!=len(outputs) or set(old)!=set(out) or len(out)!=sm['expected'] or len(out)!=sm['records']:raise ValueError('keys/count mismatch')
 for k,r in out.items():
  q=old[k];c,p=MONOMIALS[r['source']];j=r['order'];left,right=CELLS[r['cell_index']];coef=F(c*factorial(p),factorial(p-j)) if p>=j else F(0);lo=F(q['arb_lower'])+(coef*left**(p-j) if p>=j else 0);hi=F(q['arb_upper'])+(coef*right**(p-j) if p>=j else 0);pt=F(q['mp70'])
  if (F(r['mutated_lower']),F(r['mutated_upper']),F(r['point']))!=(lo,hi,pt) or r['forcing_derivative_nonzero']!=(coef!=0) or r['detected']!=(not lo<=pt<=hi):raise ValueError('arithmetic mismatch')
 if sm['detected']!=sum(x['detected'] for x in outputs) or sm['nonzero']!=sum(x['forcing_derivative_nonzero'] for x in outputs):raise ValueError('summary mismatch')
 print(json.dumps({'ok':True,'records':len(out),'nonzero':sm['nonzero'],'detected':sm['detected'],'scope':'readback arithmetic, not proof of numerical inputs'}));return 0
if __name__=='__main__':
 try:sys.exit(main())
 except (ValueError,KeyError,TypeError,OSError) as e:print(json.dumps({'ok':False,'error':str(e)}));sys.exit(1)
