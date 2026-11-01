"""Analytic stress tests: constant, sine frequencies, narrow bumps.
Synthetic error functions with exact known max 1e-4, not PINN models.
"""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,os,traceback
from flint import arb
from arb_certification import finite,exact_fraction
from arb_reconstruction import norm_enclosure
from arb_problem_families import SmoothProblem
from arb_cubic_reconstruction import SineSix,BumpSix,CubicVerifier
from arb_shaped_remainder import ShapedVerifier
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261101/analytic_shape_stress';OUT.mkdir(parents=True,exist_ok=False)
# Identical shared-data ablation formula as script88, without running its driver.
s=(ROOT/'scripts/88_matched_residual_ablation_20261101.py').read_text();s=s[s.index('def bounds('):s.index('def classification(')];env={'finite':finite,'exact_fraction':exact_fraction,'norm_enclosure':norm_enclosure,'arb':arb};exec(compile(s,'shared_envelope_bounds','exec'),env);bounds=env['bounds']
class ConstantResidual(SmoothProblem):
 def __init__(self):self.setup(103);self.boundary_error=(0,0);self.name='constant_exact'
 def residual(self,x):self.guard();return arb(8)/10000
 def r_second(self,x):self.guard();return arb(0)
 def r_third(self,x):self.guard();return arb(0)
 def r_fourth(self,x):self.guard();return arb(0)
def atom(path,obj):
 data=(json.dumps(obj,indent=2,allow_nan=False)+'\n').encode();tmp=path.with_suffix(path.suffix+'.tmp')
 with tmp.open('wb') as f:f.write(data);f.flush();os.fsync(f.fileno())
 os.replace(tmp,path)
cases=[('constant',None)]+[('sine',k) for k in (1,2,4,8,16,32)]+[('bump',w) for w in ('0.08','0.02','0.005')];grids=(16,64,128);exact=F(1,10000);rows=[]
config={'analytic_cases':cases,'grids':grids,'exact_maximum':str(exact),'methods':['global','spatial','signed_quadratic','signed_cubic'],'subdivision':2,'reference_use':'exact analytic maximum used only to check bounds and report diagnostic ratio; no selection or stopping','scope':'controlled synthetic error functions, not trained neural networks; no timing comparison','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};atom(OUT/'config.json',config)
def factory(kind,param):return ConstantResidual() if kind=='constant' else SineSix(param) if kind=='sine' else BumpSix(param)
for kind,param in cases:
 for n in grids:
  try:
   p=factory(kind,param);v=ShapedVerifier(p,'moment');v.uniform(n);q=v.certify();cnt=(v.residual_evaluations,v.curvature_evaluations);glob,local=bounds(v.prepare(v.nodes),v.depth);same=cnt==(v.residual_evaluations,v.curvature_evaluations);ql,qu=exact_fraction(q['lower']),exact_fraction(q['upper']);p3=factory(kind,param);w=CubicVerifier(p3,'moment');w.uniform(n);c=w.certify();cl,cu=exact_fraction(c['lower']),exact_fraction(c['upper'])
   for method,L,U,samples in [('global',F(0),glob,cnt[0]),('spatial',F(0),local,cnt[0]),('signed_quadratic',ql,qu,cnt[0]),('signed_cubic',cl,cu,w.residual_evaluations)]:rows.append({'case_type':kind,'parameter':param,'intervals':n,'method':method,'lower':str(L),'upper':str(U),'true_maximum':str(exact),'covers_exact':L<=exact<=U,'gap_1_02':L>0 and U<=F(51,50)*L,'upper_over_exact':float(U/exact),'residual_samples':samples,'shared_quadratic_data_no_extra':same if method!='signed_cubic' else None,'status':'finite'})
  except Exception as e:
   rows.append({'case_type':kind,'parameter':param,'intervals':n,'status':'exception','message':str(e),'traceback':traceback.format_exc()})
  atom(OUT/'records.json',rows)
 print(kind,param,'finished',flush=True)
summary={'analytic_cases':len(cases),'grids_per_case':len(grids),'expected_finite_bounds':len(cases)*len(grids)*4,'finite_bounds':sum(r['status']=='finite' for r in rows),'covers_exact':sum(r.get('covers_exact') is True for r in rows),'false_bounds':[r for r in rows if r.get('covers_exact') is False],'exceptions':[r for r in rows if r['status']=='exception'],'shared_data_no_extra':sum(r.get('shared_quadratic_data_no_extra') is True for r in rows)};atom(OUT/'summary.json',summary);print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k not in ('false_bounds','exceptions')}),'false_bounds',len(summary['false_bounds']),'exceptions',len(summary['exceptions']),flush=True)
if summary['false_bounds'] or summary['exceptions'] or summary['finite_bounds']!=summary['expected_finite_bounds']:raise SystemExit(1)
