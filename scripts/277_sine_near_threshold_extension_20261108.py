"""Post-hoc 128/256-cell follow-up on all 12 previously studied sine PINNs."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib,json,os,tempfile
from arb_certification import exact_fraction
from arb_cubic_reconstruction import SixthOrderDeep,CubicVerifier
from arb_shaped_remainder import ShapedVerifier
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results';OUT=B/'2026-11-08/sine_near_threshold_extension'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(l,u,t):return 'within' if u<=t else 'exceeds' if l>t else 'unknown'
def main():
 if OUT.exists():raise FileExistsError(OUT)
 c=json.loads((ROOT/'configs/midterm_sine_decisions_2026-11-07.json').read_text());fixture=ROOT/c['fixture'];ref=ROOT/c['reference'];summary=json.loads((ROOT/c['certificate_summary']).read_text());old=(ROOT/c['certificate_summary']).parent/'records.json'
 if sha(fixture)!=c['fixture_sha256'] or sha(ref)!=summary['references_sha256'] or sha(old)!=summary['records_sha256']:raise ValueError('input digest')
 nets=json.loads(fixture.read_text())['networks'];refs={x['network']:x for x in json.loads(ref.read_text())};oldidx={(x['network'],x['degree'],x['intervals']):x for x in json.loads(old.read_text())}
 if len(nets)!=12 or {n['label'] for n in nets}!=set(refs):raise ValueError('manifest')
 sources={str(p.relative_to(ROOT)):sha(p) for p in (ROOT/'scripts/arb_deep_network.py',ROOT/'scripts/arb_cubic_reconstruction.py',ROOT/'scripts/arb_shaped_remainder.py',ROOT/'scripts/arb_reconstruction.py',ROOT/'scripts/arb_certification.py',Path(__file__))};rows=[]
 for n,d,m in product(nets,(2,3),(128,256)):
  q=refs[n['label']];rl,ru=F(q['lower']),F(q['upper']);row={'network':n['label'],'frequency':n['frequency'],'degree':d,'mesh':m,'reference_lower':str(rl),'reference_upper':str(ru)}
  try:
   p=SixthOrderDeep(fixture,n['label'],c['fixture_sha256'],103);v=ShapedVerifier(p,'moment',subdivision=2) if d==2 else CubicVerifier(p,'moment',subdivision=2);v.uniform(m);z=v.certify();L,U=map(exact_fraction,(z['lower'],z['upper']));oldrow=oldidx[(n['label'],d,m)];mid=(rl+ru)/2
   row.update(status='complete',lower=str(L),upper=str(U),reference_contained=L<=rl<=ru<=U,matches_old=L==F(oldrow['lower']) and U==F(oldrow['upper']),decisions={str(f):dec(L,U,mid*f) for f in (F(99,100),F(101,100))})
  except Exception as e:row.update(status='error',error=type(e).__name__+': '+str(e))
  rows.append(row);print(n['label'],d,m,row['status'],flush=True)
 OUT.parent.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix='.sine_extend_',dir=OUT.parent));(stage/'records.json').write_text(json.dumps(rows,indent=2)+'\n');s={'expected':48,'records':len(rows),'complete':sum(r['status']=='complete' for r in rows),'contained':sum(r.get('reference_contained',False) for r in rows),'matches_old':sum(r.get('matches_old',False) for r in rows),'fixture_sha256':sha(fixture),'reference_sha256':sha(ref),'old_records_sha256':sha(old),'source_hashes_before_run':sources,'records_sha256':sha(stage/'records.json'),'scope':'post-hoc 1% thresholds and 128/256 meshes, 12 existing sine networks, no timing'};(stage/'summary.json').write_text(json.dumps(s,indent=2)+'\n');os.rename(stage,OUT);print(json.dumps(s));return int(s['complete']!=48 or s['contained']!=48)
if __name__=='__main__':raise SystemExit(main())
