"""Recompute selected frozen PINN certificates at three precisions."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib,json,os,tempfile
from arb_certification import exact_fraction
from arb_cubic_reconstruction import SixthOrderDeep,CubicVerifier
from arb_shaped_remainder import ShapedVerifier
ROOT=Path(__file__).resolve().parents[1];CFG=ROOT/'configs/midterm_precision_check_2026-11-08.json';BASE=ROOT/'results';OUT=BASE/'2026-11-08/midterm_precision_check'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def decision(l,u,t):return 'within' if u<=t else 'exceeds' if l>t else 'unknown'
def main():
 if OUT.exists():raise FileExistsError(OUT)
 c=json.loads(CFG.read_text());poly_script=ROOT/'scripts/74_polynomial_forcing_prototype.py';code=poly_script.read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
 if code.count(needle)!=1:raise ValueError('prototype')
 ns={'__file__':str(poly_script)};exec(compile(code.replace(needle,'OUT=None'),str(poly_script),'exec'),ns)
 old={};cases={}
 for family in ('poly','sine'):
  d=BASE/f'2026-11-07/midterm_{family}_decisions';summary=json.loads((d/'summary.json').read_text());p=d/'records.json';audit=json.loads((BASE/f'2026-11-07/midterm_{family}_decisions_readback.json').read_text())
  if sha(p)!=summary['records_sha256'] or sha(p)!=audit['records_sha256']:raise ValueError('prior digest')
  old[family]={(r['network'],r['degree'],r['mesh']):r for r in json.loads(p.read_text())}
 source={str(p.relative_to(ROOT)):sha(p) for p in (CFG,poly_script,ROOT/'scripts/arb_cubic_reconstruction.py',ROOT/'scripts/arb_shaped_remainder.py',ROOT/'scripts/arb_deep_network.py',ROOT/'scripts/arb_reconstruction.py',ROOT/'scripts/arb_certification.py',Path(__file__))}
 inputs=json.loads((ROOT/'configs/independent_poly_split_2026-11-05.json').read_text())
 for f in inputs['fixtures']:
  p=ROOT/f['path']
  if sha(p)!=f['sha256']:raise ValueError('fixture')
  for n in json.loads(p.read_text())['networks']:
   if n['label'] in c['poly_networks']:cases[('poly',n['label'])]=(n,p,f['sha256'])
 sc=json.loads((ROOT/'configs/midterm_sine_decisions_2026-11-07.json').read_text());sf=ROOT/sc['fixture']
 if sha(sf)!=sc['fixture_sha256']:raise ValueError('sine fixture')
 for n in json.loads(sf.read_text())['networks']:
  if n['label'] in c['sine_networks']:cases[('sine',n['label'])]=(n,sf,sc['fixture_sha256'])
 expected={(f,n) for f in ('poly','sine') for n in c[f+'_networks']}
 if set(cases)!=expected:raise ValueError('input manifest')
 rows=[]
 for (family,name),bits,d,m in product(sorted(cases),c['bits'],c['degrees'],c['meshes']):
  n,p,h=cases[(family,name)];prev=old[family][(name,d,m)];rl,ru=F(prev['reference_lower']),F(prev['reference_upper']);tol=F(c[family+'_tolerance']);r={'family':family,'network':name,'fixture_sha256':h,'bits':bits,'degree':d,'mesh':m,'reference_lower':str(rl),'reference_upper':str(ru),'tolerance':str(tol)}
  try:
   problem=ns['PolynomialForcing'](p,name,h,n['source_coefficients'],bits) if family=='poly' else SixthOrderDeep(p,name,h,bits)
   v=ShapedVerifier(problem,'moment',subdivision=c['bernstein_depth']) if d==2 else CubicVerifier(problem,'moment',subdivision=c['bernstein_depth']);v.uniform(m);z=v.certify();L,U=map(exact_fraction,(z['lower'],z['upper']))
   r.update(status='complete',lower=str(L),upper=str(U),reference_contained=L<=rl<=ru<=U,decision=decision(L,U,tol),reference_decision=decision(rl,ru,tol),matches_prior_103=(L==F(prev['lower']) and U==F(prev['upper'])) if bits==103 else None)
  except Exception as e:r.update(status='error',error=type(e).__name__+': '+str(e))
  rows.append(r);print(family,name,bits,d,m,r['status'],flush=True)
 OUT.parent.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix='.precision_check_',dir=OUT.parent));(stage/'locked_config.json').write_bytes(CFG.read_bytes());(stage/'records.json').write_text(json.dumps(rows,indent=2)+'\n');s={'expected':96,'records':len(rows),'complete':sum(r['status']=='complete' for r in rows),'contained':sum(r.get('reference_contained',False) for r in rows),'matches_prior_103':sum(r.get('matches_prior_103') is True for r in rows),'source_hashes_before_run':source,'previous_records_sha256':{k:sha(BASE/f'2026-11-07/midterm_{k}_decisions/records.json') for k in old},'records_sha256':sha(stage/'records.json'),'scope':c['scope']};(stage/'summary.json').write_text(json.dumps(s,indent=2)+'\n');os.rename(stage,OUT);print(json.dumps(s));return int(s['records']!=96 or s['complete']!=96 or s['contained']!=96)
if __name__=='__main__':raise SystemExit(main())
