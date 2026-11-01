"""Predefined shifted/multipeak manufactured errors; reference after runs."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from statistics import median
import json,hashlib,time,os,traceback
from mesh_problem_families import GaussianMixture,MPGaussianReference
from arb_cubic_reconstruction import CubicVerifier
from arb_shaped_remainder import ShapedVerifier
from arb_certification import exact_fraction
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261101/shifted_peak_comparison';OUT.mkdir(parents=True,exist_ok=False)
cases={'left_narrow':[{'center':'1/7','width':'1/200','weight':'1'}],'right_mid':[{'center':'5/8','width':'1/100','weight':'1'}],'near_boundary':[{'center':'49/50','width':'1/250','weight':'1'}],'two_peaks':[{'center':'1/4','width':'3/250','weight':'1/2'},{'center':'18/25','width':'1/125','weight':'1/2'}]};grid=(16,32,64,128,256,512);rows=[]
def atom(path,obj):
 d=(json.dumps(obj,indent=2,allow_nan=False)+'\n').encode();tmp=path.with_suffix(path.suffix+'.tmp')
 with tmp.open('wb') as f:f.write(d);f.flush();os.fsync(f.fileno())
 os.replace(tmp,path)
config={'cases':cases,'degrees':[2,3],'methods':['uniform','adaptive'],'budgets':grid,'repeats':3,'goal':'U/L<=1.02; cap512','reference_use':'computed only after all verification runs; no source for method selection or stopping','scope':'manufactured shifted Gaussian/mixed errors, not neural networks','timing_scope':'in process from problem construction including all trials/mesh construction; excludes imports/reference/file output','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};atom(OUT/'config.json',config)
for name,terms in cases.items():
 for deg in (2,3):
  for rep in range(1,4):
   for mesh in (('uniform','adaptive') if rep%2 else ('adaptive','uniform')):
    start=time.perf_counter();trace=[];status='budget_exhausted';error=None;nodes=None
    try:
     p=GaussianMixture(terms);v=ShapedVerifier(p,'moment') if deg==2 else CubicVerifier(p,'moment');v.uniform(grid[0])
     for n in grid:
      if mesh=='uniform':v.uniform(n)
      elif len(v.nodes)-1<n:v.refine_to(n)
      result=v.certify();L,U=exact_fraction(result['lower']),exact_fraction(result['upper']);met=L>0 and U<=F(51,50)*L;trace.append({'intervals':n,'lower':str(L),'upper':str(U),'goal_met':met,'residual_evaluations':v.residual_evaluations,'derivative_interval_evaluations':v.curvature_evaluations});nodes=[str(x) for x in v.nodes]
      if met:status='goal_met';break
    except Exception as e:status='exception';error={'message':str(e),'traceback':traceback.format_exc()}
    rows.append({'case':name,'degree':deg,'mesh':mesh,'repeat':rep,'status':status,'final':trace[-1] if trace else None,'trace':trace,'nodes':nodes,'elapsed_seconds':time.perf_counter()-start,'error':error});atom(OUT/'records.json',rows)
  print(name,'degree',deg,'done',flush=True)
refs={}
for name,terms in cases.items():
 try:refs[name]=MPGaussianReference(terms).maximum(relative=F(1,10**7),absolute=F(1,10**14),cap=5000)
 except Exception as e:refs[name]={'converged':False,'error':str(e),'traceback':traceback.format_exc()}
 atom(OUT/'references.json',refs);print('REFERENCE',name,refs[name].get('converged'),flush=True)
def fr(d):return F(int(d['numerator']),int(d['denominator']))
checks=[]
for r in rows:
 ref=refs[r['case']]
 for t in r['trace']:
  rl,ru=fr(ref['lower']),fr(ref['upper']);L,U=F(t['lower']),F(t['upper']);checks.append({'case':r['case'],'degree':r['degree'],'mesh':r['mesh'],'repeat':r['repeat'],'intervals':t['intervals'],'reference_converged':ref['converged'],'overlap':L<=ru and U>=rl,'contains_reference':L<=rl and U>=ru})
atom(OUT/'checks.json',checks);groups=[]
for name in cases:
 for deg in (2,3):
  for mesh in ('uniform','adaptive'):
   z=[r for r in rows if r['case']==name and r['degree']==deg and r['mesh']==mesh];groups.append({'case':name,'degree':deg,'mesh':mesh,'outcomes':dict(Counter(r['status'] for r in z)),'final_intervals':[r['final']['intervals'] if r['final'] else None for r in z],'median_seconds':median(r['elapsed_seconds'] for r in z)})
summary={'cases':len(cases),'requests':len(rows),'expected':len(cases)*2*2*3,'outcomes':dict(Counter(r['status'] for r in rows)),'reference_converged':sum(x['converged'] for x in refs.values()),'reference_count':len(refs),'intermediate_checks':len(checks),'reference_overlap':sum(x['overlap'] for x in checks),'reference_contained':sum(x['contains_reference'] for x in checks),'bad_overlap':[x for x in checks if not x['overlap'] or not x['reference_converged']],'exceptions':[x for x in rows if x['status']=='exception'],'groups':groups};atom(OUT/'summary.json',summary);print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k not in ('bad_overlap','exceptions','groups')}),'exceptions',len(summary['exceptions']),'bad',len(summary['bad_overlap']),flush=True)
if summary['requests']!=summary['expected'] or summary['exceptions'] or summary['bad_overlap']:raise SystemExit(1)
