"""Exact coefficient check of the polynomial forcing solution used by 291."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results/2026-11-08';OUT=B/'polynomial_solution_exact_check.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def derivative(p):return [i*p[i] for i in range(1,len(p))]
def main():
 if OUT.exists():raise FileExistsError(OUT)
 config=json.loads((ROOT/'configs/midterm_independent_sample_2026-11-08.json').read_text());fixture=ROOT/config['training_dir']/'frozen_networks.json'
 if sha(fixture)!=config['fixture_sha256']:raise ValueError('fixture hash')
 nets=json.loads(fixture.read_text())['networks'];rows=[]
 for n in nets:
  coeff=list(map(F,n['source_coefficients']));u=[F(0)]*(len(coeff)+2)
  for j,c in enumerate(coeff):
   d=F((j+1)*(j+2));u[1]+=c/d;u[j+2]-=c/d
  second=[-x for x in derivative(derivative(u))];zero_left=u[0]==0;zero_right=sum(u)==0
  if second!=coeff or not zero_left or not zero_right:raise ValueError('forcing/boundary identity')
  rows.append({'network':n['label'],'source':n['source_name'],'forcing_coefficients':[str(x) for x in coeff],'solution_coefficients':[str(x) for x in u],'forcing_exact':True,'zero_boundaries_exact':True})
 if len(rows)!=12 or len({r['network'] for r in rows})!=12:raise ValueError('manifest')
 result={'networks':12,'forcing_and_boundaries_exact':12,'fixture_sha256':sha(fixture),'script_sha256':sha(Path(__file__)),'rows':rows,'scope':'exact rational polynomial identities only; does not certify neural-network arithmetic'}
 with OUT.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps({k:v for k,v in result.items() if k!='rows'}))
if __name__=='__main__':main()
