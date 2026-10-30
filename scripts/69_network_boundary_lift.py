"""Known affine boundary perturbations of frozen networks; reference uses known sine solution."""
from pathlib import Path
from fractions import Fraction as F
import importlib.util,json,hashlib,time
from arb_cubic_reconstruction import SixthOrderDeep,CubicVerifier
from arb_shaped_remainder import ShapedVerifier
from arb_certification import exact_fraction
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261030/network_boundary_lift';OUT.mkdir(parents=True,exist_ok=False)
source=(ROOT/'scripts/67_boundary_lift_prototype.py').read_text().split('boundary=')[0];source=source.replace("OUT=ROOT.parent/'run_20261030/boundary_lift';OUT.mkdir(parents=True,exist_ok=False)","OUT=None");ns={'__file__':str(ROOT/'scripts/67_boundary_lift_prototype.py')};exec(compile(source,'lift_prototype','exec'),ns);lift_norm=ns['lift_norm']
spec=importlib.util.spec_from_file_location('mpref24',ROOT/'scripts/24_depth_continuous_reference.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class LiftedReference(m.MPDeepReference):
 def __init__(self,row,g0,g1):super().__init__(row);self.g0=g0;self.g1=g1
 def error(self,x):return super().error(x)+self.point(self.g0)+self.point(self.g1-self.g0)*self.point(x)
def frac(d):return F(int(d['numerator']),int(d['denominator']))
fixture=ROOT/'configs/frozen_networks_holdout_2026-10-30.json';data=fixture.read_bytes();sha=hashlib.sha256(data).hexdigest();rows=[];refs={};cases=['depth1_k2_seed49','depth3_k5_seed49'];boundaries={'zero':(F(0),F(0)),'small_opposite':(F(-1,10000),F(1,10000)),'asymmetric':(F(1,1000),F(-1,2000))}
for name in cases:
 raw=next(x for x in json.loads(data)['networks'] if x['label']==name)
 for bname,(g0,g1) in boundaries.items():
  ref=LiftedReference(raw,g0,g1).maximum(relative=F(1,10**7),absolute=F(1,10**13),cap=5000);refs[name+'__'+bname]=ref;(OUT/'references.json').write_text(json.dumps(refs,indent=2)+'\n')
  rl,ru=frac(ref['lower']),frac(ref['upper']);print('REFERENCE',name,bname,ref['converged'],ref['splits'],flush=True)
  for degree in (2,3):
   for count in (128,256):
    prob=SixthOrderDeep(fixture,name,sha,103);v=ShapedVerifier(prob,'moment') if degree==2 else CubicVerifier(prob,'moment');nodes=[F(i,count) for i in range(count+1)];cells=v.prepare(nodes);zl,zu=lift_norm(v,cells,g0,g1);_,_,du=v._norms(cells);lo=exact_fraction(max(ns['arb'](0),(zl-du).lower()));hi=exact_fraction((zu+du).upper());contained=lo<=rl and hi>=ru;overlap=lo<=ru and hi>=rl
    row={'network':name,'boundary':bname,'g0':str(g0),'g1':str(g1),'degree':degree,'intervals':count,'lower':str(lo),'upper':str(hi),'reference_lower':str(rl),'reference_upper':str(ru),'reference_converged':ref['converged'],'reference_overlap':overlap,'reference_contained':contained};rows.append(row);(OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n');print(name,bname,degree,count,'contains',contained,flush=True)
summary={'networks':len(cases),'reference_cases':len(refs),'converged':sum(x['converged'] for x in refs.values()),'certificates':len(rows),'overlap':sum(x['reference_overlap'] for x in rows),'contained':sum(x['reference_contained'] for x in rows),'failures':[x for x in rows if not x['reference_overlap'] or not x['reference_converged']],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'reference_only':'known sine solution and mpmath interval backend; certifier uses neither'};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k!='failures'}),'failures',len(summary['failures']),flush=True)
