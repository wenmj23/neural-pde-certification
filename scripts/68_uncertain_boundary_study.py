"""Known boundary intervals: affine midpoint certificate plus endpoint radius.
Research-only synthetic tests; residual must remain certified on the domain.
"""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261030/uncertain_boundary';OUT.mkdir(parents=True,exist_ok=False)
source=(ROOT/'scripts/67_boundary_lift_prototype.py').read_text().split('boundary=')[0];source=source.replace("OUT=ROOT.parent/'run_20261030/boundary_lift';OUT.mkdir(parents=True,exist_ok=False)","OUT=None");ns={'__file__':str(ROOT/'scripts/67_boundary_lift_prototype.py')};exec(compile(source,'lift_prototype','exec'),ns)
P=ns['POLYS'];certify=ns['certify'];reference=ns['reference'];names=['zero','constant','linear','node_hidden_quartic','asymmetric_sextic']; nominal={'zero':(F(0),F(0)),'opposite':(F(-1,10000),F(1,10000)),'unequal':(F(1,20000),F(-1,40000))};radii={'exact':(F(0),F(0)),'small':(F(1,1000000),F(1,2000000)),'large':(F(1,10000),F(1,20000))};meshes={'four':[F(i,4) for i in range(5)],'nonuniform':[F(0),F(1,17),F(1,3),F(7,10),F(1)]};rows=[];refs={}
for name in names:
 for bname,(g0,g1) in nominal.items():
  for rname,(r0,r1) in radii.items():
   key=(name,bname,rname); scenarios=[(g0-r0,g1-r1),(g0-r0,g1+r1),(g0+r0,g1-r1),(g0+r0,g1+r1),(g0,g1)]
   refs[key]=[(a,b,*reference(P[name],a,b,depth=10)) for a,b in scenarios]
   for mesh,nodes in meshes.items():
    for degree in (2,3):
     lo,hi,_,_=certify(P[name],nodes,degree,g0,g1);eta=max(r0,r1);L=max(F(0),lo-eta);U=hi+eta
     checks=[{'g0':str(a),'g1':str(b),'lower':str(rl),'upper':str(ru),'compatible':L<=ru and U>=rl} for a,b,rl,ru in refs[key]]
     rows.append({'case':name,'nominal':bname,'radius_case':rname,'radius0':str(r0),'radius1':str(r1),'mesh':mesh,'degree':degree,'nominal_L':str(lo),'nominal_U':str(hi),'uncertain_L':str(L),'uncertain_U':str(U),'scenarios':checks,'all_compatible':all(x['compatible'] for x in checks),'exact_radius_unchanged':(L==lo and U==hi) if eta==0 else None});(OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n')
 print(name,'done',flush=True)
summary={'families':len(names),'certificates':len(rows),'scenario_comparisons':sum(len(x['scenarios']) for x in rows),'scenario_passed':sum(z['compatible'] for x in rows for z in x['scenarios']),'exact_radius_checks':sum(x['exact_radius_unchanged'] is not None for x in rows),'exact_radius_passed':sum(x['exact_radius_unchanged'] is True for x in rows),'failures':[x for x in rows if not x['all_compatible'] or x['exact_radius_unchanged'] is False],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k!='failures'}),'failures',len(summary['failures']),flush=True)
