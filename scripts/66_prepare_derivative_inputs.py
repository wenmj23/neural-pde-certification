"""Create derivative-audit input descriptors from repository fixtures."""
from pathlib import Path
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
CASES=['depth1_k1_seed7','depth1_k5_seed17','depth2_k1_seed17','depth2_k5_seed27','depth3_k1_seed27','depth3_k5_seed7','depth1_k2_seed49','depth3_k2_seed49','depth3_k5_seed49']
FIXTURES=['frozen_networks_depth_2026-10-29.json','frozen_networks_selection_fresh_2026-10-29.json','frozen_networks_holdout_2026-10-30.json']
p=argparse.ArgumentParser();p.add_argument('--dest',type=Path,required=True);arg=p.parse_args();arg.dest.mkdir(parents=True,exist_ok=False);found={}
for fn in FIXTURES:
 path=ROOT/'configs'/fn; raw=path.read_bytes();sha=hashlib.sha256(raw).hexdigest()
 for n in json.loads(raw)['networks']:
  if n['label'] in CASES:
   if n['label'] in found:raise ValueError('Ambiguous network')
   found[n['label']]={'network':{'fixture':str(path.resolve()),'label':n['label'],'sha256':sha}}
if set(found)!=set(CASES):raise ValueError('Missing network')
for name,c in found.items():(arg.dest/(name+'.case.json')).write_text(json.dumps(c,indent=2)+'\n')
print('Prepared',len(found),'descriptors')
