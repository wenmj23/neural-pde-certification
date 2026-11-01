"""Same one-sided tolerance task; quadratic global/spatial/signed and cubic.
Exploratory thresholds chosen after earlier experiments, not a blind evaluation.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from statistics import median
import json,hashlib,time,os,traceback
from flint import arb
from arb_certification import finite,exact_fraction
from arb_reconstruction import norm_enclosure
from arb_shaped_remainder import ShapedVerifier
from arb_cubic_reconstruction import CubicVerifier
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261101/matched_acceptance_benchmark';OUT.mkdir(parents=True,exist_ok=False)
s=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];s=s.replace("OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)","OUT=None");ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(s,'poly_forcing_definitions','exec'),ns)
fixture=ROOT.parent/'run_20261101/poly_new_training/frozen_polynomial_pinns.json';data=fixture.read_bytes();sha=hashlib.sha256(data).hexdigest();assert sha==(fixture.parent/'fixture_sha256.txt').read_text().strip();nets=json.loads(data)['networks']
methods=('global','spatial','signed_quadratic','signed_cubic');tols=(F(1,1000),F(1,100000));grid=(16,32,64,128);rows=[]
def atomic(path,value):
 data=(json.dumps(value,indent=2,allow_nan=False)+'\n').encode();tmp=path.with_suffix(path.suffix+'.tmp')
 with tmp.open('wb') as f:f.write(data);f.flush();os.fsync(f.fileno())
 os.replace(tmp,path)

def upper_for_cells(cells,depth,method):
 env=[]
 for c in cells:
  r0,rm,r1=c['r0'],c['rm'],c['r1'];_,p=norm_enclosure([r0,2*rm-(r0+r1)/2,r1],depth);env.append(finite(p+c['delta'],'rho').upper())
 if method=='global':return exact_fraction(finite(max(env)/8,'global').upper())
 slope=arb(0)
 for c,rho in zip(cells,env):slope=finite(slope+rho*c['h']*(1-(c['a']+c['b'])/2),'slope')
 val=upper=arb(0)
 for c,rho in zip(cells,env):
  h=c['h'];coeff=[val,val+slope*h/2,val+slope*h-rho*h*h/2];_,b=norm_enclosure(coeff,depth);upper=max(upper,b);val=coeff[-1];slope=finite(slope-rho*h,'slope update')
 if not val.contains(0):raise ArithmeticError('envelope boundary mismatch')
 return exact_fraction(finite(upper,'spatial upper').upper())
config={'fixture_sha256':sha,'methods':methods,'tolerances':[str(t) for t in tols],'uniform_mesh_sequence':grid,'repetitions':3,'stop':'only if upper <= tolerance; otherwise unknown after budget','timing_scope':'in-process elapsed from network object creation through final certificate; includes all earlier meshes; excludes imports, reference, writes','reference_use':'read only after all requests; no solver stopping/selection','study_limit':'same existing networks and thresholds chosen after previous experiments; exploratory, not blind','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};atomic(OUT/'config.json',config)
for net in nets:
 for tau in tols:
  for rep in range(3):
   order=methods[rep:]+methods[:rep]
   for method in order:
    tick=time.perf_counter();trace=[];status='unknown';error=None
    try:
     p=ns['PolynomialForcing'](fixture,net['label'],sha,net['source_coefficients']);v=CubicVerifier(p,'moment') if method=='signed_cubic' else ShapedVerifier(p,'moment')
     for n in grid:
      v.uniform(n)
      if method in ('global','spatial'):U=upper_for_cells(v.prepare(v.nodes),v.depth,method)
      else:U=exact_fraction(v.certify()['upper'])
      met=U<=tau;trace.append({'intervals':n,'upper':str(U),'within':met,'residual_evaluations':v.residual_evaluations,'derivative_evaluations':v.curvature_evaluations})
      if met:status='within';break
    except Exception as e:status='exception';error={'message':str(e),'traceback':traceback.format_exc()}
    rows.append({'network':net['label'],'source':net['source_name'],'depth':net['depth'],'seed':net['seed'],'tolerance':str(tau),'repeat':rep+1,'method':method,'status':status,'final':trace[-1] if trace else None,'trace':trace,'elapsed_seconds':time.perf_counter()-tick,'error':error})
  atomic(OUT/'records.json',rows)
 print(net['label'],'finished',flush=True)
refs=json.loads((ROOT.parent/'run_20261101/poly_new_certification/references.json').read_text());fr=ns['fr'];checks=[]
for r in rows:
 ref=refs[r['network']];lo,hi=fr(ref['lower']),fr(ref['upper']);tau=F(r['tolerance']);refclass='within' if hi<=tau else 'exceeds' if lo>tau else 'reference_unknown';upper=F(r['final']['upper']) if r['final'] else None
 checks.append({'network':r['network'],'tolerance':r['tolerance'],'repeat':r['repeat'],'method':r['method'],'reference_class':refclass,'reference_converged':ref['converged'],'upper_covers_reference_lower':bool(upper is not None and upper>=lo),'false_accept':r['status']=='within' and refclass=='exceeds'})
atomic(OUT/'checks.json',checks)
groups=[]
for tau in tols:
 for method in methods:
  subset=[r for r in rows if r['tolerance']==str(tau) and r['method']==method];eligible=[r for r in subset if next(c for c in checks if (c['network'],c['tolerance'],c['repeat'],c['method'])==(r['network'],r['tolerance'],r['repeat'],r['method']))['reference_class']=='within'];groups.append({'tolerance':str(tau),'method':method,'requests':len(subset),'outcomes':dict(Counter(r['status'] for r in subset)),'eligible_requests':len(eligible),'eligible_success':sum(r['status']=='within' for r in eligible),'median_success_seconds':median(r['elapsed_seconds'] for r in eligible if r['status']=='within') if any(r['status']=='within' for r in eligible) else None})
pairs=[]
for net in nets:
 for tau in tols:
  for a,b in (('global','signed_quadratic'),('spatial','signed_quadratic'),('signed_quadratic','signed_cubic')):
   aa=[r for r in rows if r['network']==net['label'] and r['tolerance']==str(tau) and r['method']==a];bb=[r for r in rows if r['network']==net['label'] and r['tolerance']==str(tau) and r['method']==b]
   if len(aa)==len(bb)==3 and all(r['status']=='within' for r in aa+bb):
    am,bm=median(r['elapsed_seconds'] for r in aa),median(r['elapsed_seconds'] for r in bb);pairs.append({'network':net['label'],'tolerance':str(tau),'baseline':a,'comparison':b,'baseline_median_seconds':am,'comparison_median_seconds':bm,'comparison_over_baseline':bm/am})
summary={'networks':len(nets),'case_tolerance_pairs':len(nets)*len(tols),'requests':len(rows),'expected_requests':len(nets)*len(tols)*len(methods)*3,'reference_checks':len(checks),'reference_converged':sum(c['reference_converged'] for c in checks),'bounds_cover_reference_lower':sum(c['upper_covers_reference_lower'] for c in checks),'false_accepts':[c for c in checks if c['false_accept']],'exceptions':[r for r in rows if r['status']=='exception'],'groups':groups,'matched_success_pairs':pairs,'repeats_are_independent_networks':False};atomic(OUT/'summary.json',summary)
print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k not in ('groups','matched_success_pairs','exceptions','false_accepts')}),'exceptions',len(summary['exceptions']),'false_accepts',len(summary['false_accepts']),flush=True)
if summary['requests']!=summary['expected_requests'] or summary['exceptions'] or summary['false_accepts'] or summary['bounds_cover_reference_lower']!=len(checks):raise SystemExit(1)
