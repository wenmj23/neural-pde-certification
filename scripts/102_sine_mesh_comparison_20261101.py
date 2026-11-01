"""Synthetic sine errors: uniform vs residual-driven adaptive mesh."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from statistics import median
import json,hashlib,time,os,traceback
from arb_cubic_reconstruction import SineSix,CubicVerifier
from arb_shaped_remainder import ShapedVerifier
from arb_certification import exact_fraction
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261101/sine_mesh_comparison';OUT.mkdir(parents=True,exist_ok=False)
widths=(1,2,4,8,16,32);grid=(16,32,64,128,256,512);true=F(1,10000);rows=[]
def atom(path,obj):
 data=(json.dumps(obj,indent=2,allow_nan=False)+'\n').encode();tmp=path.with_suffix(path.suffix+'.tmp')
 with tmp.open('wb') as f:f.write(data);f.flush();os.fsync(f.fileno())
 os.replace(tmp,path)
config={'cases':widths,'degree':[2,3],'mesh':['uniform','adaptive'],'budgets':grid,'repeats':3,'goal':'U/L<=1.02; maximum 512 intervals','exact_maximum':str(true),'reference_use':'known exact error used only after each result to check inclusion; never used by refinement or stopping','sampling':'synthetic sine error, not neural network; same frequencies as fixed-grid screen','timing_scope':'in process from problem creation through all trials and mesh refinement, excludes import, file writing','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};atom(OUT/'config.json',config)
for width in widths:
 for deg in (2,3):
  for rep in range(1,4):
   for mesh in (('uniform','adaptive') if rep%2 else ('adaptive','uniform')):
    start=time.perf_counter();trace=[];status='budget_exhausted';error=None;nodes=None
    try:
     p=SineSix(width);v=ShapedVerifier(p,'moment') if deg==2 else CubicVerifier(p,'moment');v.uniform(grid[0])
     for n in grid:
      if mesh=='uniform':v.uniform(n)
      elif len(v.nodes)-1<n:v.refine_to(n)
      r=v.certify();L,U=exact_fraction(r['lower']),exact_fraction(r['upper']);met=L>0 and U<=F(51,50)*L;trace.append({'intervals':n,'lower':str(L),'upper':str(U),'goal_met':met,'covers_exact':L<=true<=U,'residual_evaluations':v.residual_evaluations,'derivative_interval_evaluations':v.curvature_evaluations});nodes=[str(x) for x in v.nodes]
      if met:status='goal_met';break
    except Exception as e:status='exception';error={'message':str(e),'traceback':traceback.format_exc()}
    rows.append({'frequency':width,'degree':deg,'mesh':mesh,'repeat':rep,'status':status,'final':trace[-1] if trace else None,'trace':trace,'nodes':nodes,'elapsed_seconds':time.perf_counter()-start,'error':error});atom(OUT/'records.json',rows)
  print('frequency',width,'degree',deg,'finished',flush=True)
groups=[]
for width in widths:
 for degree in (2,3):
  for mesh in ('uniform','adaptive'):
   z=[r for r in rows if r['frequency']==width and r['degree']==degree and r['mesh']==mesh];groups.append({'frequency':width,'degree':degree,'mesh':mesh,'outcomes':dict(Counter(r['status'] for r in z)),'final_intervals':[r['final']['intervals'] if r['final'] else None for r in z],'median_seconds':median(r['elapsed_seconds'] for r in z)})
checks=[t for r in rows for t in r['trace']];summary={'requests':len(rows),'expected':len(widths)*2*2*3,'outcomes':dict(Counter(r['status'] for r in rows)),'intermediate_bounds':len(checks),'exact_coverage':sum(t['covers_exact'] for t in checks),'exceptions':[r for r in rows if r['status']=='exception'],'bad_bounds':[r for r in rows if any(not t['covers_exact'] for t in r['trace'])],'groups':groups};atom(OUT/'summary.json',summary);print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k not in ('groups','exceptions','bad_bounds')}),'exceptions',len(summary['exceptions']),'bad',len(summary['bad_bounds']),flush=True)
if summary['requests']!=summary['expected'] or summary['exceptions'] or summary['bad_bounds']:raise SystemExit(1)
