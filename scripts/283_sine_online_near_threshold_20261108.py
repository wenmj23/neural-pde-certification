"""Actual sequential mesh strategy on previously studied sine candidates."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib,json,math,os,random,tempfile,time
from arb_certification import exact_fraction
from arb_cubic_reconstruction import SixthOrderDeep,CubicVerifier
from arb_shaped_remainder import ShapedVerifier
ROOT=Path(__file__).resolve().parents[1];CFG=ROOT/'configs/sine_online_near_threshold_2026-11-08.json';OUT=ROOT/'results/2026-11-08/sine_online_near_threshold'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(l,u,t):return 'within' if u<=t else 'exceeds' if l>t else 'unknown'
def main():
 if OUT.exists():raise FileExistsError(OUT)
 c=json.loads(CFG.read_text());base=ROOT/c['source_config'];sc=json.loads(base.read_text());fixture=ROOT/sc['fixture'];ref=ROOT/sc['reference'];prior_dir=(ROOT/sc['certificate_summary']).parent;prior=json.loads((prior_dir/'summary.json').read_text());old=prior_dir/'records.json'
 if sha(fixture)!=sc['fixture_sha256'] or sha(ref)!=prior['references_sha256'] or sha(old)!=prior['records_sha256']:raise ValueError('input hash')
 nets={n['label']:n for n in json.loads(fixture.read_text())['networks']};refs={q['network']:q for q in json.loads(ref.read_text())};oldidx={(x['network'],x['degree'],x['intervals']):x for x in json.loads(old.read_text())}
 if len(nets)!=12 or set(nets)!=set(refs):raise ValueError('network manifest')
 paths=[CFG,base,ROOT/'scripts/arb_deep_network.py',ROOT/'scripts/arb_cubic_reconstruction.py',ROOT/'scripts/arb_shaped_remainder.py',ROOT/'scripts/arb_reconstruction.py',ROOT/'scripts/arb_certification.py',Path(__file__)];sources={str(p.relative_to(ROOT)):sha(p) for p in paths};requests=list(product(sorted(nets),c['degrees'],c['threshold_factors'],c['strategies']));random.Random(c['shuffle_seed']).shuffle(requests);rows=[]
 for name,degree,fac,strategy in requests:
  q=refs[name];rl,ru=F(q['lower']),F(q['upper']);tol=(rl+ru)*F(fac)/2;attempts=[];row={'network':name,'degree':degree,'factor':fac,'strategy':strategy,'reference_lower':str(rl),'reference_upper':str(ru),'tolerance':str(tol)}
  try:
   start=time.perf_counter();p=SixthOrderDeep(fixture,name,sc['fixture_sha256'],c['precision_bits']);v=ShapedVerifier(p,'moment',subdivision=c['bernstein_depth']) if degree==2 else CubicVerifier(p,'moment',subdivision=c['bernstein_depth']);core=0.0
   for mesh in c['strategies'][strategy]:
    v.uniform(mesh);tick=time.perf_counter();z=v.certify();elapsed=time.perf_counter()-tick;core+=elapsed;L,U=map(exact_fraction,(z['lower'],z['upper']));o=oldidx[(name,degree,mesh)];value=dec(L,U,tol)
    attempts.append({'mesh':mesh,'seconds':elapsed,'lower':str(L),'upper':str(U),'reference_contained':L<=rl<=ru<=U,'matches_old':L==F(o['lower']) and U==F(o['upper']),'decision':value,'residual_calls_cumulative':v.residual_evaluations,'derivative_calls_cumulative':v.curvature_evaluations})
    if value!='unknown':break
   inclusive=time.perf_counter()-start
   if not all(math.isfinite(z) and z>0 for z in (core,inclusive)) or inclusive<core:raise ValueError('timing')
   row.update(status='complete',decision=attempts[-1]['decision'],attempts=attempts,core_seconds=core,inclusive_seconds=inclusive)
  except Exception as e:row.update(status='error',error=type(e).__name__+': '+str(e),attempts=attempts)
  rows.append(row);print(name,degree,fac,strategy,row['status'],flush=True)
 OUT.parent.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix='.online_near_',dir=OUT.parent));(stage/'locked_config.json').write_bytes(CFG.read_bytes());(stage/'records.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n');s={'expected':len(requests),'records':len(rows),'complete':sum(r['status']=='complete' for r in rows),'attempts':sum(len(r['attempts']) for r in rows),'contained':sum(a['reference_contained'] for r in rows for a in r['attempts']),'matches_old':sum(a['matches_old'] for r in rows for a in r['attempts']),'inputs_sha256':{'fixture':sha(fixture),'reference':sha(ref),'old_records':sha(old)},'source_hashes_before_run':sources,'records_sha256':sha(stage/'records.json'),'scope':c['scope'],'measurement':c['measurement']};(stage/'summary.json').write_text(json.dumps(s,indent=2)+'\n');os.rename(stage,OUT);print(json.dumps(s));return int(s['records']!=96 or s['complete']!=96 or s['contained']!=s['attempts'])
if __name__=='__main__':raise SystemExit(main())
