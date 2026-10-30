"""Point-reference diagnostics for bounded uncertain network boundary errors."""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib
from flint import arb
from arb_cubic_reconstruction import SixthOrderDeep,CubicVerifier
from arb_certification import exact_fraction
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261030/network_uncertain_boundary';OUT.mkdir(parents=True,exist_ok=False)
source=(ROOT/'scripts/69_network_boundary_lift.py').read_text().split("fixture=ROOT/'configs/")[0];source=source.replace("OUT=ROOT.parent/'run_20261030/network_boundary_lift';OUT.mkdir(parents=True,exist_ok=False)","OUT=None");ns={'__file__':str(ROOT/'scripts/69_network_boundary_lift.py')};exec(compile(source,'network_lift_definitions','exec'),ns);lift_norm=ns['lift_norm'];Reference=ns['LiftedReference'];frac=ns['frac']
fixture=ROOT/'configs/frozen_networks_holdout_2026-10-30.json';data=fixture.read_bytes();sha=hashlib.sha256(data).hexdigest();allnet=json.loads(data)['networks'];names=['depth1_k2_seed49','depth3_k5_seed49'];g0=g1=F(1,1000);radii={'small':(F(1,100000),F(1,200000)),'large':(F(1,2000),F(1,4000))};rows=[]; refs={}
def decide(l,u,t):return 'within' if u<=t else 'exceeds' if l>t else 'unknown'
for name in names:
 raw=next(x for x in allnet if x['label']==name)
 p=SixthOrderDeep(fixture,name,sha,103);v=CubicVerifier(p,'moment');nodes=[F(i,256) for i in range(257)];cells=v.prepare(nodes);zl,zu=lift_norm(v,cells,g0,g1);_,_,du=v._norms(cells);l=exact_fraction(max(arb(0),(zl-du).lower()));u=exact_fraction((zu+du).upper())
 for label,(r0,r1) in radii.items():
  eta=max(r0,r1);L=max(F(0),l-eta);U=u+eta;scenarios=[(g0-r0,g1-r1),(g0-r0,g1+r1),(g0+r0,g1-r1),(g0+r0,g1+r1),(g0,g1)];checks=[]
  for i,(a,b) in enumerate(scenarios):
   ref=Reference(raw,a,b).maximum(relative=F(1,10**7),absolute=F(1,10**13),cap=5000);key=f'{name}__{label}__{i}';refs[key]=ref;(OUT/'references.json').write_text(json.dumps(refs,indent=2)+'\n');rl,ru=frac(ref['lower']),frac(ref['upper']);checks.append({'g0':str(a),'g1':str(b),'reference_lower':str(rl),'reference_upper':str(ru),'converged':ref['converged'],'contained':L<=rl and U>=ru,'overlap':L<=ru and U>=rl});print(name,label,i,'contained',checks[-1]['contained'],flush=True)
  dec={str(t):{'center':decide(l,u,t),'uncertain':decide(L,U,t)} for t in (F(1,1000),F(1,100))};rows.append({'network':name,'radius':label,'g0_center':str(g0),'g1_center':str(g1),'radius0':str(r0),'radius1':str(r1),'center_lower':str(l),'center_upper':str(u),'uncertain_lower':str(L),'uncertain_upper':str(U),'decisions':dec,'checks':checks});(OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n')
summary={'networks':len(names),'uncertain_certificates':len(rows),'scenario_references':len(refs),'converged':sum(x['converged'] for x in refs.values()),'contained':sum(c['contained'] for x in rows for c in x['checks']),'overlap':sum(c['overlap'] for x in rows for c in x['checks']),'lost_decisions':sum(v['center']!='unknown' and v['uncertain']=='unknown' for x in rows for v in x['decisions'].values()),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps(summary),flush=True)
