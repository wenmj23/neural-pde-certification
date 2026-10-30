"""Fixed public-network portfolio; compare public CLI degrees, no training."""
import hashlib, json, platform, subprocess, sys, time
from fractions import Fraction as F
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT.parent/'run_20261030'/'frozen_portfolio';OUT.mkdir(parents=True,exist_ok=False)
fixtures=['frozen_networks_depth_2026-10-29.json','frozen_networks_selection_fresh_2026-10-29.json','frozen_networks_holdout_2026-10-30.json']
rows=[]; refs=json.loads((ROOT.parent/'run_20261030/public_crosscheck/published_references.json').read_text())
def fr(d):return F(int(d['numerator']),int(d['denominator']))
for fixture in fixtures:
 path=ROOT/'configs'/fixture; sha=hashlib.sha256(path.read_bytes()).hexdigest(); data=json.loads(path.read_text())
 for network in data['networks']:
  name=network['label'];k=network['frequency']; in_path=OUT/(name+'.case.json')
  case={'schema_version':1,'name':name,'problem':{'operator':'-u_dd','source':'(k*pi)^2*sin(k*pi*x)','frequency':k,'domain':[0,1],'boundary':[0,0]},'network':{'fixture':str(path),'label':name,'sha256':sha,'output_transform':'x*(1-x)*raw_network(x)'},'verification':{'degree':3,'mesh':'uniform','goal':'enclosure','gap_ratio':'1.02','error_tolerance':'0.001','initial_intervals':16,'max_intervals':1024,'max_attempts':8,'time_limit_seconds':'20','precision_bits':103}}
  in_path.write_text(json.dumps(case,indent=2)+'\n'); pair=[]
  for degree in (2,3):
   dst=OUT/f'{name}__degree{degree}'; cmd=[sys.executable,str(ROOT/'scripts/certify.py'),'--case',str(in_path),'--out',str(dst),'--degree',str(degree)]
   t=time.perf_counter()
   try:
    proc=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,timeout=40); ret=proc.returncode; stdout=proc.stdout;stderr=proc.stderr
   except subprocess.TimeoutExpired as e:
    ret='external_timeout';stdout=str(e.stdout);stderr=str(e.stderr)
   p=dst/'result.json';d=json.loads(p.read_text()) if p.exists() else None;enc=d['enclosure'] if d else None; ref=refs.get(name)
   ref_ok=(fr(enc['lower'])<=fr(ref['lower']) and fr(ref['upper'])<=fr(enc['upper'])) if enc and ref else None
   row={'name':name,'fixture':fixture,'depth':network['depth'],'frequency':k,'degree':degree,'returncode':ret,'stop_reason':d['stop_reason'] if d else 'no_result','goal_met':d['goal_met'] if d else False,'decision':d['accuracy_decision'] if d else None,'intervals':enc['intervals'] if enc else None,'lower':str(fr(enc['lower'])) if enc else None,'upper':str(fr(enc['upper'])) if enc else None,'reference_contained':ref_ok,'wall_external_seconds':time.perf_counter()-t,'stdout':stdout,'stderr':stderr}
   rows.append(row);pair.append(row);(OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n');print(name,degree,row['stop_reason'],row['intervals'],flush=True)
  if all(x['lower'] is not None for x in pair):
   overlap=max(F(x['lower']) for x in pair)<=min(F(x['upper']) for x in pair)
  else:overlap=None
  (OUT/(name+'.overlap.json')).write_text(json.dumps({'name':name,'overlap':overlap},indent=2)+'\n')
summary={'commit':'1948ea20ae6d88395dc417c97581a78a427484da','networks':len(rows)//2,'requests':len(rows),'goal_met':sum(x['goal_met'] for x in rows),'stop_reasons':dict(Counter(x['stop_reason'] for x in rows)),'unexpected_exit':[r for r in rows if r['returncode'] not in (0,2)],'reference_checked':sum(x['reference_contained'] is not None for x in rows),'reference_contained':sum(x['reference_contained'] is True for x in rows),'overlaps':[json.loads(p.read_text()) for p in OUT.glob('*.overlap.json')],'python':platform.python_version(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n'); print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k not in ('overlaps','unexpected_exit')}),'nonoverlap',sum(x['overlap'] is False for x in summary['overlaps']),'unexpected',len(summary['unexpected_exit']),flush=True)
