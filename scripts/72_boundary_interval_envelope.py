"""Research-only corner upper and boundary endpoint lower for uncertain errors.
Uses convexity of the infinity norm in affine boundary parameters and certified
bounds for each corner; residual is fixed as the boundary lift varies.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import json,hashlib
from flint import arb
from arb_cubic_reconstruction import SixthOrderDeep,CubicVerifier
from arb_certification import exact_fraction
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261030/boundary_interval_envelope';OUT.mkdir(parents=True,exist_ok=False)
source=(ROOT/'scripts/67_boundary_lift_prototype.py').read_text().split('boundary=')[0];source=source.replace("OUT=ROOT.parent/'run_20261030/boundary_lift';OUT.mkdir(parents=True,exist_ok=False)","OUT=None");ns={'__file__':str(ROOT/'scripts/67_boundary_lift_prototype.py')};exec(compile(source,'lift_definition','exec'),ns);lift=ns['lift_norm']
fixture=ROOT/'configs/frozen_networks_holdout_2026-10-30.json';data=fixture.read_bytes();sha=hashlib.sha256(data).hexdigest();nets=json.loads(data)['networks'];references=json.loads((ROOT.parent/'run_20261030/network_uncertain_boundary/references.json').read_text())
BOXES={'equal_small':(F(1,1000),F(1,1000),F(1,100000),F(1,200000)), 'equal_large':(F(1,1000),F(1,1000),F(1,2000),F(1,4000)), 'opposite_large':(F(-1,1000),F(2,1000),F(1,2000),F(1,4000)), 'positive_large':(F(5,1000),F(5,1000),F(2,1000),F(1,1000))}
def classify(l,u,t):return 'within' if u<=t else 'exceeds' if l>t else 'unknown'
def frac(d):return F(int(d['numerator']),int(d['denominator']))
def min_abs(a,b):return F(0) if a<=0<=b else min(abs(a),abs(b))
rows=[]
for net in nets:
 name=net['label'];problem=SixthOrderDeep(fixture,name,sha,103);v=CubicVerifier(problem,'moment');nodes=[F(i,256) for i in range(257)];cells=v.prepare(nodes);_,_,du=v._norms(cells)
 def certified(g0,g1):
  zl,zu=lift(v,cells,g0,g1);return exact_fraction(max(arb(0),(zl-du).lower())),exact_fraction((zu+du).upper())
 for box,(c0,c1,r0,r1) in BOXES.items():
  center=certified(c0,c1);eta=max(r0,r1);oldL=max(F(0),center[0]-eta);oldU=center[1]+eta
  corners=[(c0-r0,c1-r1),(c0-r0,c1+r1),(c0+r0,c1-r1),(c0+r0,c1+r1)];cs=[certified(a,b) for a,b in corners];cornerU=max(u for l,u in cs);endL=max(min_abs(c0-r0,c0+r0),min_abs(c1-r1,c1+r1));newL=max(oldL,endL);newU=min(oldU,cornerU)
  comparisons=[]
  if name in ('depth1_k2_seed49','depth3_k5_seed49') and box in ('equal_small','equal_large'):
   lab='small' if box=='equal_small' else 'large'
   for i in range(5):
    ref=references[f'{name}__{lab}__{i}'];rl,ru=frac(ref['lower']),frac(ref['upper']);comparisons.append({'reference_key':f'{name}__{lab}__{i}','converged':ref['converged'],'contained':newL<=rl and ru<=newU,'overlap':newL<=ru and rl<=newU})
  decisions=[]
  for tau in (F(1,1000),F(1,100)):
   old=classify(oldL,oldU,tau);new=classify(newL,newU,tau);decisions.append({'tolerance':str(tau),'old':old,'new':new,'contradiction':old!='unknown' and new!='unknown' and old!=new})
  row={'network':name,'box':box,'center':[str(c0),str(c1)],'radius':[str(r0),str(r1)],'center_L':str(center[0]),'center_U':str(center[1]),'baseline_L':str(oldL),'baseline_U':str(oldU),'endpoint_L':str(endL),'corner_U':str(cornerU),'refined_L':str(newL),'refined_U':str(newU),'corners':[{'g0':str(a),'g1':str(b),'lower':str(l),'upper':str(u)} for (a,b),(l,u) in zip(corners,cs)],'valid_order':newL<=newU,'subset_baseline':oldL<=newL<=newU<=oldU,'improved_lower':newL>oldL,'improved_upper':newU<oldU,'decisions':decisions,'reference_checks':comparisons}
  rows.append(row);(OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n');print(name,box,'lower',row['improved_lower'],'upper',row['improved_upper'],'decisions',[(z['old'],z['new']) for z in decisions],flush=True)
refs=[c for r in rows for c in r['reference_checks']];dec=[d for r in rows for d in r['decisions']];summary={'networks':len(nets),'boxes_per_network':len(BOXES),'certificates':len(rows),'valid_order':sum(r['valid_order'] for r in rows),'subset_baseline':sum(r['subset_baseline'] for r in rows),'improved_lower':sum(r['improved_lower'] for r in rows),'improved_upper':sum(r['improved_upper'] for r in rows),'reference_checks':len(refs),'reference_contained':sum(x['contained'] for x in refs),'reference_overlaps':sum(x['overlap'] for x in refs),'decision_pairs':len(dec),'baseline_decided':sum(x['old']!='unknown' for x in dec),'refined_decided':sum(x['new']!='unknown' for x in dec),'gained_decision':sum(x['old']=='unknown' and x['new']!='unknown' for x in dec),'contradictions':sum(x['contradiction'] for x in dec),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps(summary),flush=True)
