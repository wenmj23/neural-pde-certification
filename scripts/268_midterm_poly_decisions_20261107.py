"""Fresh certification of old networks under fixed meshes and thresholds."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib,json,os,tempfile
from arb_certification import exact_fraction
from arb_cubic_reconstruction import CubicVerifier
from arb_shaped_remainder import ShapedVerifier
ROOT=Path(__file__).resolve().parents[1];CFG=ROOT/'configs/midterm_poly_decisions_2026-11-07.json';OUT=ROOT/'results/2026-11-07/midterm_poly_decisions'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def refval(d):return F(int(d['numerator']),int(d['denominator']))
def dec(l,u,t):return 'within' if u<=t else 'exceeds' if l>t else 'unknown'
def main():
 if OUT.exists():raise FileExistsError(OUT)
 c=json.loads(CFG.read_text());inputs=ROOT/c['fixtures_config'];input_data=json.loads(inputs.read_text());poly=ROOT/'scripts/74_polynomial_forcing_prototype.py';code=poly.read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
 if code.count(needle)!=1:raise ValueError('prototype')
 ns={'__file__':str(poly)};exec(compile(code.replace(needle,'OUT=None'),str(poly),'exec'),ns)
 provenance={str(p.relative_to(ROOT)):sha(p) for p in (CFG,inputs,poly,ROOT/'scripts/arb_cubic_reconstruction.py',ROOT/'scripts/arb_shaped_remainder.py',ROOT/'scripts/arb_reconstruction.py',ROOT/'scripts/arb_certification.py',Path(__file__))}
 cases=[]
 for f in input_data['fixtures']:
  path,ref=ROOT/f['path'],ROOT/f['reference']
  if sha(path)!=f['sha256'] or sha(ref)!=f['reference_sha256']:raise ValueError('fixture/reference hash')
  refs=json.loads(ref.read_text());nets=json.loads(path.read_text())['networks']
  if len(nets)!=12 or {n['label'] for n in nets}!=set(refs):raise ValueError('manifest')
  for n in nets:
   q=refs[n['label']]
   if not q['converged']:raise ValueError('reference')
   cases.append((n,path,f['sha256'],refval(q['lower']),refval(q['upper'])))
 if len(cases)!=24 or len({n['label'] for n,_,_,_,_ in cases})!=24:raise ValueError('duplicate')
 rows=[]
 for (n,path,h,rl,ru),d,m in product(cases,c['degrees'],c['meshes']):
  row={'network':n['label'],'source':n['source_name'],'fixture_sha256':h,'degree':d,'mesh':m,'reference_lower':str(rl),'reference_upper':str(ru)}
  try:
   p=ns['PolynomialForcing'](path,n['label'],h,n['source_coefficients'],c['precision_bits']);v=ShapedVerifier(p,'moment',subdivision=c['bernstein_depth']) if d==2 else CubicVerifier(p,'moment',subdivision=c['bernstein_depth']);v.uniform(m);r=v.certify();L,U=map(exact_fraction,(r['lower'],r['upper']))
   row.update(status='complete',lower=str(L),upper=str(U),reference_contained=L<=rl<=ru<=U,decisions={t:dec(L,U,F(t)) for t in c['tolerances']},reference_decisions={t:dec(rl,ru,F(t)) for t in c['tolerances']})
  except Exception as e:row.update(status='error',error=type(e).__name__+': '+str(e))
  rows.append(row);print(n['label'],d,m,row['status'],flush=True)
 OUT.parent.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix='.midterm_decisions_',dir=OUT.parent));(stage/'locked_config.json').write_bytes(CFG.read_bytes());(stage/'records.json').write_text(json.dumps(rows,indent=2)+'\n');s={'expected':144,'records':len(rows),'complete':sum(r['status']=='complete' for r in rows),'reference_contained':sum(r.get('reference_contained',False) for r in rows),'source_hashes_before_run':provenance,'records_sha256':sha(stage/'records.json'),'scope':c['scope']};(stage/'summary.json').write_text(json.dumps(s,indent=2)+'\n');os.rename(stage,OUT);print(json.dumps(s));return int(s['records']!=144 or s['complete']!=144 or s['reference_contained']!=144)
if __name__=='__main__':raise SystemExit(main())
