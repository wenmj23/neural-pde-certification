"""Fixed-case comparison: affine reconstruction vs triangle bound, research prototype."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import json,hashlib
from flint import arb
from arb_cubic_reconstruction import SixthOrderDeep,CubicVerifier
from arb_certification import exact_fraction
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261030/boundary_network_comparison';OUT.mkdir(parents=True,exist_ok=False)
source=(ROOT/'scripts/69_network_boundary_lift.py').read_text().split("fixture=ROOT/'configs/")[0];source=source.replace("OUT=ROOT.parent/'run_20261030/network_boundary_lift';OUT.mkdir(parents=True,exist_ok=False)","OUT=None");ns={'__file__':str(ROOT/'scripts/69_network_boundary_lift.py')};exec(compile(source,'network_lift_definitions','exec'),ns);lift_norm=ns['lift_norm'];Reference=ns['LiftedReference'];frac=ns['frac']
fixture=ROOT/'configs/frozen_networks_holdout_2026-10-30.json';data=fixture.read_bytes();sha=hashlib.sha256(data).hexdigest();networks=json.loads(data)['networks']; boundaries={'none':(F(0),F(0)),'small_opposite':(F(-1,10000),F(1,10000)),'equal_positive':(F(1,1000),F(1,1000)),'large_asymmetric':(F(-5,1000),F(2,1000))};tol=[F(1,1000),F(1,100)];rows=[];refs={}
def classify(l,u,t):return 'within' if u<=t else 'exceeds' if l>t else 'unknown'
for raw in networks:
 name=raw['label']
 for bname,(g0,g1) in boundaries.items():
  ref=Reference(raw,g0,g1).maximum(relative=F(1,10**7),absolute=F(1,10**13),cap=5000);refs[name+'__'+bname]=ref;(OUT/'references.json').write_text(json.dumps(refs,indent=2)+'\n');rl,ru=frac(ref['lower']),frac(ref['upper']);print('REFERENCE',name,bname,'converged',ref['converged'],flush=True)
  for n in (128,256):
   problem=SixthOrderDeep(fixture,name,sha,103);v=CubicVerifier(problem,'moment');nodes=[F(i,n) for i in range(n+1)];cells=v.prepare(nodes);zl,zu=lift_norm(v,cells,g0,g1);_,_,du=v._norms(cells);L=exact_fraction(max(arb(0),(zl-du).lower()));U=exact_fraction((zu+du).upper());base=v.certify(nodes);bl,bu=exact_fraction(base['lower']),exact_fraction(base['upper']);A=max(abs(g0),abs(g1));TL=max(F(0),bl-A,A-bu);TU=bu+A
   decisions=[]
   for t in tol:
    d=classify(L,U,t);b=classify(TL,TU,t);decisions.append({'tolerance':str(t),'direct':d,'triangle':b,'opposed':d!='unknown' and b!='unknown' and d!=b})
   row={'network':name,'boundary':bname,'g0':str(g0),'g1':str(g1),'intervals':n,'direct_lower':str(L),'direct_upper':str(U),'triangle_lower':str(TL),'triangle_upper':str(TU),'reference_lower':str(rl),'reference_upper':str(ru),'reference_converged':ref['converged'],'direct_reference_contained':L<=rl and U>=ru,'triangle_reference_contained':TL<=rl and TU>=ru,'direct_reference_overlap':L<=ru and U>=rl,'triangle_reference_overlap':TL<=ru and TU>=rl,'ranges_overlap':max(L,TL)<=min(U,TU),'width_relation':'narrower' if U-L<TU-TL else 'equal' if U-L==TU-TL else 'wider','decisions':decisions}
   rows.append(row);(OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n');print(name,bname,n,row['width_relation'],[z['direct']+'/'+z['triangle'] for z in decisions],flush=True)
dec=[d for r in rows for d in r['decisions']];summary={'network_count':len(networks),'boundary_conditions':len(boundaries),'reference_cases':len(refs),'reference_converged':sum(r['converged'] for r in refs.values()),'certificates':len(rows),'direct_reference_contained':sum(r['direct_reference_contained'] for r in rows),'triangle_reference_contained':sum(r['triangle_reference_contained'] for r in rows),'direct_overlap':sum(r['direct_reference_overlap'] for r in rows),'triangle_overlap':sum(r['triangle_reference_overlap'] for r in rows),'range_overlap':sum(r['ranges_overlap'] for r in rows),'width_relation':dict(Counter(r['width_relation'] for r in rows)),'decision_pairs':len(dec),'direct_determined':sum(x['direct']!='unknown' for x in dec),'triangle_determined':sum(x['triangle']!='unknown' for x in dec),'direct_only':sum(x['direct']!='unknown' and x['triangle']=='unknown' for x in dec),'triangle_only':sum(x['triangle']!='unknown' and x['direct']=='unknown' for x in dec),'opposed':sum(x['opposed'] for x in dec),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps(summary),flush=True)
