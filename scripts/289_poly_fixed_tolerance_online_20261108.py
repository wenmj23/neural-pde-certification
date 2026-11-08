"""Repeat online early stopping for two fixed tolerances on polynomial networks."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib,json,math,os,random,tempfile,time
from arb_certification import exact_fraction
from arb_cubic_reconstruction import CubicVerifier
from arb_shaped_remainder import ShapedVerifier
ROOT=Path(__file__).resolve().parents[1];CFG=ROOT/'configs/poly_fixed_tolerance_online_2026-11-08.json';B=ROOT/'results/2026-11-08';OUT=B/'poly_fixed_tolerance_online'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def val(x):return F(int(x['numerator']),int(x['denominator']))
def dec(l,u,t):return 'within' if u<=t else 'exceeds' if l>t else 'unknown'
def main():
 if OUT.exists():raise FileExistsError(OUT)
 c=json.loads(CFG.read_text());base=ROOT/c['source_config'];bc=json.loads(base.read_text());ic=ROOT/bc['fixtures_config'];inputs=json.loads(ic.read_text());prior=ROOT/'results/2026-11-07/midterm_poly_decisions/records.json';audit=json.loads((ROOT/'results/2026-11-07/midterm_poly_decisions_readback.json').read_text())
 if sha(prior)!=audit['records_sha256']:raise ValueError('prior hash')
 old={(r['network'],r['degree'],r['mesh']):r for r in json.loads(prior.read_text())};poly=ROOT/'scripts/74_polynomial_forcing_prototype.py';code=poly.read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
 if code.count(needle)!=1:raise ValueError('prototype')
 ns={'__file__':str(poly)};exec(compile(code.replace(needle,'OUT=None'),str(poly),'exec'),ns);nets={}
 for f in inputs['fixtures']:
  path,ref=ROOT/f['path'],ROOT/f['reference']
  if sha(path)!=f['sha256'] or sha(ref)!=f['reference_sha256']:raise ValueError('input hash')
  refs=json.loads(ref.read_text())
  for n in json.loads(path.read_text())['networks']:
   q=refs[n['label']]
   if n['label'] in nets or not q['converged']:raise ValueError('network/ref')
   nets[n['label']]=(n,path,f['sha256'],val(q['lower']),val(q['upper']))
 if len(nets)!=24:raise ValueError('manifest')
 files=[CFG,base,ic,poly,ROOT/'scripts/arb_deep_network.py',ROOT/'scripts/arb_cubic_reconstruction.py',ROOT/'scripts/arb_shaped_remainder.py',ROOT/'scripts/arb_reconstruction.py',ROOT/'scripts/arb_certification.py',Path(__file__)];sources={str(p.relative_to(ROOT)):sha(p) for p in files};requests=list(product(sorted(nets),c['degrees'],c['tolerances'],c['strategies']));random.Random(c['shuffle_seed']).shuffle(requests);rows=[]
 for name,degree,t,strategy in requests:
  n,path,h,rl,ru=nets[name];tol=F(t);attempts=[];row={'network':name,'source':n['source_name'],'fixture_sha256':h,'degree':degree,'tolerance':t,'strategy':strategy,'reference_lower':str(rl),'reference_upper':str(ru),'reference_decision':dec(rl,ru,tol)}
  try:
   start=time.perf_counter();p=ns['PolynomialForcing'](path,name,h,n['source_coefficients'],c['precision_bits']);v=ShapedVerifier(p,'moment',subdivision=c['bernstein_depth']) if degree==2 else CubicVerifier(p,'moment',subdivision=c['bernstein_depth']);core=0.0
   for mesh in c['strategies'][strategy]:
    v.uniform(mesh);tick=time.perf_counter();z=v.certify();elapsed=time.perf_counter()-tick;core+=elapsed;L,U=map(exact_fraction,(z['lower'],z['upper']));o=old[(name,degree,mesh)];value=dec(L,U,tol)
    attempts.append({'mesh':mesh,'seconds':elapsed,'lower':str(L),'upper':str(U),'reference_contained':L<=rl<=ru<=U,'matches_old':L==F(o['lower']) and U==F(o['upper']),'decision':value,'residual_calls_cumulative':v.residual_evaluations,'derivative_calls_cumulative':v.curvature_evaluations})
    if value!='unknown':break
   inclusive=time.perf_counter()-start
   if not all(math.isfinite(x) and x>0 for x in (core,inclusive)) or inclusive<core:raise ValueError('timing')
   row.update(status='complete',decision=attempts[-1]['decision'],attempts=attempts,core_seconds=core,inclusive_seconds=inclusive)
  except Exception as e:row.update(status='error',error=type(e).__name__+': '+str(e),attempts=attempts)
  rows.append(row);print(name,degree,t,strategy,row['status'],flush=True)
 B.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix='.poly_fixed_',dir=B));(stage/'locked_config.json').write_bytes(CFG.read_bytes());(stage/'records.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n');s={'expected':len(requests),'records':len(rows),'complete':sum(r['status']=='complete' for r in rows),'attempts':sum(len(r['attempts']) for r in rows),'contained':sum(a['reference_contained'] for r in rows for a in r['attempts']),'matches_old':sum(a['matches_old'] for r in rows for a in r['attempts']),'prior_records_sha256':sha(prior),'source_hashes_before_run':sources,'records_sha256':sha(stage/'records.json'),'scope':c['scope'],'measurement':c['measurement']};(stage/'summary.json').write_text(json.dumps(s,indent=2)+'\n');os.rename(stage,OUT);print(json.dumps(s));return int(s['records']!=192 or s['complete']!=192 or s['contained']!=s['attempts'])
if __name__=='__main__':raise SystemExit(main())
