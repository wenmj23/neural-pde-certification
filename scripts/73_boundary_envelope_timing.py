"""Paired wall-time diagnostics; fixed 256 mesh, includes all common setup."""
from pathlib import Path
from fractions import Fraction as F
from statistics import median
import json,hashlib,time
from flint import arb
from arb_cubic_reconstruction import SixthOrderDeep,CubicVerifier
from arb_certification import exact_fraction
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261030/boundary_envelope_timing';OUT.mkdir(parents=True,exist_ok=False)
source=(ROOT/'scripts/67_boundary_lift_prototype.py').read_text().split('boundary=')[0];source=source.replace("OUT=ROOT.parent/'run_20261030/boundary_lift';OUT.mkdir(parents=True,exist_ok=False)","OUT=None");ns={'__file__':str(ROOT/'scripts/67_boundary_lift_prototype.py')};exec(compile(source,'lift_definition','exec'),ns);lift=ns['lift_norm']
fixture=ROOT/'configs/frozen_networks_holdout_2026-10-30.json';data=fixture.read_bytes();sha=hashlib.sha256(data).hexdigest();nets=json.loads(data)['networks'];prev=json.loads((ROOT.parent/'run_20261030/boundary_interval_envelope/records.json').read_text());boxes={'equal_small':(F(1,1000),F(1,1000),F(1,100000),F(1,200000)),'equal_large':(F(1,1000),F(1,1000),F(1,2000),F(1,4000))};rows=[]
def min_abs(a,b):return F(0) if a<=0<=b else min(abs(a),abs(b))
def run(name,box,method):
 c0,c1,r0,r1=boxes[box];start=time.perf_counter();problem=SixthOrderDeep(fixture,name,sha,103);v=CubicVerifier(problem,'moment');nodes=[F(i,256) for i in range(257)];cells=v.prepare(nodes);_,_,du=v._norms(cells)
 def bound(g0,g1):
  zl,zu=lift(v,cells,g0,g1);return exact_fraction(max(arb(0),(zl-du).lower())),exact_fraction((zu+du).upper())
 center=bound(c0,c1);eta=max(r0,r1);L=max(F(0),center[0]-eta);U=center[1]+eta
 if method=='refined':
  corner=[bound(a,b) for a in (c0-r0,c0+r0) for b in (c1-r1,c1+r1)];L=max(L,min_abs(c0-r0,c0+r0),min_abs(c1-r1,c1+r1));U=min(U,max(x[1] for x in corner))
 return L,U,time.perf_counter()-start
for n in nets:
 name=n['label']
 for box in boxes:
  prior=next(r for r in prev if r['network']==name and r['box']==box)
  for rep in range(3):
   for method in (('baseline','refined') if rep%2==0 else ('refined','baseline')):
    lo,hi,elapsed=run(name,box,method);expectedL=F(prior['baseline_L' if method=='baseline' else 'refined_L']);expectedU=F(prior['baseline_U' if method=='baseline' else 'refined_U']);ok=(lo==expectedL and hi==expectedU)
    rows.append({'network':name,'box':box,'repeat':rep+1,'method':method,'elapsed_seconds':elapsed,'exact_match_prior':ok});(OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n')
  print(name,box,'done',flush=True)
groups=[]
for n in nets:
 for box in boxes:
  q=[r for r in rows if r['network']==n['label'] and r['box']==box];base=median(r['elapsed_seconds'] for r in q if r['method']=='baseline');ref=median(r['elapsed_seconds'] for r in q if r['method']=='refined');prior=next(r for r in prev if r['network']==n['label'] and r['box']==box);old=F(prior['baseline_U']);new=F(prior['refined_U']);groups.append({'network':n['label'],'box':box,'baseline_median_s':base,'refined_median_s':ref,'time_ratio':ref/base,'upper_reduction_percent':float((old-new)/old*100)})
summary={'networks':len(nets),'boxes':len(boxes),'paired_configurations':len(groups),'trials':len(rows),'exact_match_prior':sum(r['exact_match_prior'] for r in rows),'median_of_median_time_ratios':median(x['time_ratio'] for x in groups),'min_time_ratio':min(x['time_ratio'] for x in groups),'max_time_ratio':max(x['time_ratio'] for x in groups),'groups':groups,'timing_scope':'fresh fixture+problem setup, 256-cell prepare, D and all lift calculations; excludes Python import, file writing, reference; three paired repeats on one machine','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k!='groups'}),flush=True)
