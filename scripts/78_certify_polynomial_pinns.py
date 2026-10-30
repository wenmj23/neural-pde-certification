"""Post-training polynomial PINN study. Reference used only for validation."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import hashlib,json
from arb_cubic_reconstruction import CubicVerifier
from arb_shaped_remainder import ShapedVerifier
from arb_certification import exact_fraction
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261030/polynomial_pinn_certificates';OUT.mkdir(parents=True,exist_ok=False)
s=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];s=s.replace("OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)","OUT=None");ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(s,'poly_forcing_definitions','exec'),ns)
fixture=ROOT.parent/'run_20261030/polynomial_pinn_training/frozen_polynomial_pinns.json';data=fixture.read_bytes();sha=hashlib.sha256(data).hexdigest();assert sha==(fixture.parent/'fixture_sha256.txt').read_text().strip();nets=json.loads(data)['networks'];rows=[];refs={};fr=ns['fr']
for n in nets:
 name=n['label'];coeff=n['source_coefficients'];ref=ns['PolynomialReference'](n,coeff).maximum(relative=F(1,10**7),absolute=F(1,10**14),cap=5000);refs[name]=ref;(OUT/'references.json').write_text(json.dumps(refs,indent=2)+'\n');rl,ru=fr(ref['lower']),fr(ref['upper']);print('REFERENCE',name,'converged',ref['converged'],'splits',ref['splits'],flush=True)
 for degree in (2,3):
  for count in (32,64,128,256):
   p=ns['PolynomialForcing'](fixture,name,sha,coeff);v=ShapedVerifier(p,'moment') if degree==2 else CubicVerifier(p,'moment');v.uniform(count);r=v.certify();l,u=exact_fraction(r['lower']),exact_fraction(r['upper']);row={'network':name,'source':n['source_name'],'seed':n['seed'],'degree':degree,'intervals':count,'lower':str(l),'upper':str(u),'remainder_upper':str(exact_fraction(r['remainder_upper'])),'reference_lower':str(rl),'reference_upper':str(ru),'reference_converged':ref['converged'],'reference_overlap':l<=ru and u>=rl,'reference_contained':l<=rl and u>=ru,'gap_1_02':l>0 and u<=F(51,50)*l,'decision_0_001':'within' if u<=F(1,1000) else 'exceeds' if l>F(1,1000) else 'unknown'};rows.append(row);(OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n');print(name,degree,count,'tight',row['gap_1_02'],'contains',row['reference_contained'],flush=True)
summary={'fixture_sha256':sha,'networks':len(nets),'sources':len(set(x['source'] for x in rows)),'reference_cases':len(refs),'reference_converged':sum(x['converged'] for x in refs.values()),'certificates':len(rows),'reference_overlap':sum(x['reference_overlap'] for x in rows),'reference_contained':sum(x['reference_contained'] for x in rows),'gap_1_02':sum(x['gap_1_02'] for x in rows),'decisions':dict(Counter(x['decision_0_001'] for x in rows)),'failures':[x for x in rows if not x['reference_overlap'] or not x['reference_converged']],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'reference_role':'validation after freeze only'};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k!='failures'}),'failures',len(summary['failures']),flush=True)
