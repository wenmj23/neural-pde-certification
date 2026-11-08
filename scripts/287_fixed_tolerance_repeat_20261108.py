"""Repeated randomized online timing for one fixed absolute tolerance."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib,json,math,os,random,tempfile,time
from arb_certification import exact_fraction
from arb_cubic_reconstruction import SixthOrderDeep,CubicVerifier
from arb_shaped_remainder import ShapedVerifier
ROOT=Path(__file__).resolve().parents[1];CFG=ROOT/'configs/sine_fixed_tolerance_repeat_2026-11-08.json';B=ROOT/'results/2026-11-08';OUT=B/'sine_fixed_tolerance_repeat'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(l,u,t):return 'within' if u<=t else 'exceeds' if l>t else 'unknown'
def main():
 if OUT.exists():raise FileExistsError(OUT)
 c=json.loads(CFG.read_text());base=ROOT/c['base_config'];bc=json.loads(base.read_text());sc=json.loads((ROOT/bc['source_config']).read_text());fixture=ROOT/sc['fixture'];ref=ROOT/sc['reference'];old_dir=(ROOT/sc['certificate_summary']).parent;old=old_dir/'records.json';orig=B/'sine_fixed_tolerance_online/records.json';audit=json.loads((B/'sine_fixed_tolerance_online_readback.json').read_text());oldsum=json.loads((old_dir/'summary.json').read_text())
 if sha(orig)!=audit['records_sha256'] or sha(fixture)!=sc['fixture_sha256'] or sha(ref)!=oldsum['references_sha256'] or sha(old)!=oldsum['records_sha256']:raise ValueError('input hash')
 nets={n['label']:n for n in json.loads(fixture.read_text())['networks']};refs={q['network']:q for q in json.loads(ref.read_text())};origidx={(r['network'],r['degree'],r['tolerance'],r['strategy']):r for r in json.loads(orig.read_text())};oldidx={(r['network'],r['degree'],r['intervals']):r for r in json.loads(old.read_text())}
 if len(nets)!=12 or set(nets)!=set(refs):raise ValueError('network manifest')
 paths=[CFG,base,ROOT/bc['source_config'],ROOT/'scripts/arb_deep_network.py',ROOT/'scripts/arb_cubic_reconstruction.py',ROOT/'scripts/arb_shaped_remainder.py',ROOT/'scripts/arb_reconstruction.py',ROOT/'scripts/arb_certification.py',Path(__file__)];sources={str(p.relative_to(ROOT)):sha(p) for p in paths};combos=list(product(sorted(nets),bc['degrees'],bc['strategies']));rng=random.Random(c['shuffle_seed']);rows=[];tol=F(c['tolerance'])
 for phase,repeats in (('warmup',c['warmups']),('timed',c['timed_repeats'])):
  for iteration in range(repeats):
   order=combos[:];rng.shuffle(order)
   for name,degree,strategy in order:
    q=refs[name];rl,ru=F(q['lower']),F(q['upper']);attempts=[];r={'network':name,'degree':degree,'tolerance':str(tol),'strategy':strategy,'phase':phase,'iteration':iteration,'reference_lower':str(rl),'reference_upper':str(ru)}
    try:
     start=time.perf_counter();p=SixthOrderDeep(fixture,name,sc['fixture_sha256'],bc['precision_bits']);v=ShapedVerifier(p,'moment',subdivision=bc['bernstein_depth']) if degree==2 else CubicVerifier(p,'moment',subdivision=bc['bernstein_depth']);core=0.0
     for mesh in bc['strategies'][strategy]:
      v.uniform(mesh);tick=time.perf_counter();z=v.certify();dt=time.perf_counter()-tick;core+=dt;L,U=map(exact_fraction,(z['lower'],z['upper']));value=dec(L,U,tol);o=oldidx[(name,degree,mesh)]
      attempts.append({'mesh':mesh,'seconds':dt,'lower':str(L),'upper':str(U),'decision':value,'reference_contained':L<=rl<=ru<=U,'matches_old':L==F(o['lower']) and U==F(o['upper'])})
      if value!='unknown':break
     total=time.perf_counter()-start;prior=origidx[(name,degree,str(tol),strategy)];same=attempts[-1]['decision']==prior['decision'] and [a['mesh'] for a in attempts]==[a['mesh'] for a in prior['attempts']]
     if not all(math.isfinite(x) and x>0 for x in (core,total)) or total<core:raise ValueError('timing')
     r.update(status='complete',decision=attempts[-1]['decision'],attempts=attempts,core_seconds=core,inclusive_seconds=total,path_matches_previous=same)
    except Exception as e:r.update(status='error',error=type(e).__name__+': '+str(e),attempts=attempts)
    rows.append(r)
   print(phase,iteration,'completed',flush=True)
 B.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix='.fixed_repeat_',dir=B));(stage/'locked_config.json').write_bytes(CFG.read_bytes());(stage/'records.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n');attempt_count=sum(len(r['attempts']) for r in rows);s={'expected':len(combos)*(c['warmups']+c['timed_repeats']),'records':len(rows),'complete':sum(r['status']=='complete' for r in rows),'attempts':attempt_count,'contained':sum(a['reference_contained'] for r in rows for a in r['attempts']),'matches_old':sum(a['matches_old'] for r in rows for a in r['attempts']),'paths_match_previous':sum(r.get('path_matches_previous',False) for r in rows),'inputs_sha256':{'fixture':sha(fixture),'reference':sha(ref),'old_records':sha(old),'online_records':sha(orig)},'source_hashes_before_run':sources,'records_sha256':sha(stage/'records.json'),'scope':c['scope'],'measurement':c['measurement']};(stage/'summary.json').write_text(json.dumps(s,indent=2)+'\n');os.rename(stage,OUT);print(json.dumps(s));return int(s['records']!=192 or s['complete']!=192 or s['contained']!=attempt_count or s['paths_match_previous']!=192)
if __name__=='__main__':raise SystemExit(main())
