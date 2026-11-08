"""Fresh fixed-budget sine PINN decisions against saved numerical references."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib,json,os,tempfile
from arb_certification import exact_fraction
from arb_cubic_reconstruction import SixthOrderDeep,CubicVerifier
from arb_shaped_remainder import ShapedVerifier
ROOT=Path(__file__).resolve().parents[1];CFG=ROOT/'configs/midterm_sine_decisions_2026-11-07.json';OUT=ROOT/'results/2026-11-07/midterm_sine_decisions'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(l,u,t):return 'within' if u<=t else 'exceeds' if l>t else 'unknown'
def main():
 if OUT.exists():raise FileExistsError(OUT)
 c=json.loads(CFG.read_text());fixture=ROOT/c['fixture'];ref=ROOT/c['reference'];prior=json.loads((ROOT/c['certificate_summary']).read_text())
 if sha(fixture)!=c['fixture_sha256'] or sha(ref)!=prior['references_sha256'] or prior['fixture_sha256']!=c['fixture_sha256']:raise ValueError('input digest')
 nets=json.loads(fixture.read_text())['networks'];refs={r['network']:r for r in json.loads(ref.read_text())}
 if len(nets)!=12 or len(refs)!=12 or {n['label'] for n in nets}!=set(refs):raise ValueError('manifest')
 sources={str(p.relative_to(ROOT)):sha(p) for p in (CFG,ROOT/'scripts/arb_deep_network.py',ROOT/'scripts/arb_cubic_reconstruction.py',ROOT/'scripts/arb_shaped_remainder.py',ROOT/'scripts/arb_reconstruction.py',ROOT/'scripts/arb_certification.py',Path(__file__))};rows=[]
 for n,d,m in product(nets,c['degrees'],c['meshes']):
  q=refs[n['label']];rl,ru=F(q['lower']),F(q['upper']);row={'network':n['label'],'frequency':n['frequency'],'depth':n['depth'],'seed':n['seed'],'degree':d,'mesh':m,'reference_lower':str(rl),'reference_upper':str(ru)}
  try:
   if not q['converged']:raise ValueError('reference unconverged')
   p=SixthOrderDeep(fixture,n['label'],c['fixture_sha256'],c['precision_bits']);v=ShapedVerifier(p,'moment',subdivision=c['bernstein_depth']) if d==2 else CubicVerifier(p,'moment',subdivision=c['bernstein_depth']);v.uniform(m);r=v.certify();L,U=map(exact_fraction,(r['lower'],r['upper']))
   row.update(status='complete',lower=str(L),upper=str(U),reference_contained=L<=rl<=ru<=U,decisions={t:dec(L,U,F(t)) for t in c['tolerances']},reference_decisions={t:dec(rl,ru,F(t)) for t in c['tolerances']})
  except Exception as e:row.update(status='error',error=type(e).__name__+': '+str(e))
  rows.append(row);print(n['label'],d,m,row['status'],flush=True)
 OUT.parent.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix='.sine_decisions_',dir=OUT.parent));(stage/'locked_config.json').write_bytes(CFG.read_bytes());(stage/'records.json').write_text(json.dumps(rows,indent=2)+'\n');s={'expected':72,'records':len(rows),'complete':sum(r['status']=='complete' for r in rows),'reference_contained':sum(r.get('reference_contained',False) for r in rows),'reference_sha256':sha(ref),'source_hashes_before_run':sources,'records_sha256':sha(stage/'records.json'),'scope':c['scope']};(stage/'summary.json').write_text(json.dumps(s,indent=2)+'\n');os.rename(stage,OUT);print(json.dumps(s));return int(s['records']!=72 or s['complete']!=72 or s['reference_contained']!=72)
if __name__=='__main__':raise SystemExit(main())
