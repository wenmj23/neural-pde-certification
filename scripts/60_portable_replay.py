"""Replay saved CLI snapshots from an unrelated working directory."""
import json,subprocess,sys,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; BASE=ROOT.parent/'run_20261030/precision_crosscheck';OUT=ROOT.parent/'run_20261030/portable_replay';OUT.mkdir(parents=True,exist_ok=False);empty=OUT/'empty_cwd';empty.mkdir();rows=[]
names=['depth1_k1_seed7','depth1_k5_seed37','depth2_k1_seed17','depth3_k5_seed7','depth1_k2_seed49','depth3_k5_seed49']
for name in names:
 src=BASE/f'{name}__p103';old=json.loads((src/'result.json').read_text());dest=OUT/name;cmd=[sys.executable,str(src/'source_snapshot/certify.py'),'--case',str(src/'case.json'),'--out',str(dest)]
 p=subprocess.run(cmd,cwd=empty,capture_output=True,text=True,timeout=60); result=dest/'result.json';new=json.loads(result.read_text()) if result.exists() else None
 keys=('lower','upper','z_lower','z_upper','remainder_upper');equal=bool(new and all(old['enclosure'][k]['numerator']==new['enclosure'][k]['numerator'] and old['enclosure'][k]['denominator']==new['enclosure'][k]['denominator'] for k in keys)); snapshot=bool(new and old['source_sha256']==new['source_sha256']);ok=bool(p.returncode==0 and equal and snapshot and old['accuracy_decision']==new['accuracy_decision'])
 row={'name':name,'returncode':p.returncode,'exact_endpoints_equal':equal,'source_hashes_equal':snapshot,'passed':ok,'stdout':p.stdout,'stderr':p.stderr};rows.append(row);(OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n');print(name,'passed',ok,flush=True)
summary={'replays':len(rows),'passed':sum(r['passed'] for r in rows),'failed':[r for r in rows if not r['passed']],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',summary,flush=True)
