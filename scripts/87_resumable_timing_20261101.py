"""Paired full-process time comparison: original writer vs journal prototype."""
from pathlib import Path
import hashlib,json,os,shutil,statistics,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261101/resumable_timing';OUT.mkdir(parents=True,exist_ok=False);runs=[]
base=json.loads((ROOT.parent/'run_20261101/poly_fixed_tolerance/records.json').read_text());strip=lambda x:{k:v for k,v in x.items() if k!='elapsed_seconds'}
for trial in range(1,4):
 work=OUT/f'work_{trial}';repo=work/'repo';repo.mkdir(parents=True);shutil.copytree(ROOT/'scripts',repo/'scripts',ignore=shutil.ignore_patterns('__pycache__'))
 inp=work/'run_20261101/poly_new_training';inp.mkdir(parents=True);src=ROOT.parent/'run_20261101/poly_new_training'
 for n in ('frozen_polynomial_pinns.json','fixture_sha256.txt'):shutil.copy2(src/n,inp/n)
 refdir=work/'run_20261101/poly_new_certification';refdir.mkdir();shutil.copy2(ROOT.parent/'run_20261101/poly_new_certification/references.json',refdir/'references.json')
 for method in (('legacy','journal') if trial%2 else ('journal','legacy')):
  out=work/'run_20261101/poly_fixed_tolerance' if method=='legacy' else work/'journal_output'
  cmd=[sys.executable,str(repo/'scripts/82_fixed_tolerance_poly_20261101.py')] if method=='legacy' else [sys.executable,str(repo/'scripts/85_resumable_decisions_20261101.py'),'--out',str(out)]
  start=time.perf_counter();p=subprocess.run(cmd,cwd=repo,capture_output=True,text=True,timeout=180);elapsed=time.perf_counter()-start;f=out/'records.json';r=json.loads(f.read_text()) if f.exists() else [];equal=len(r)==len(base)==72 and all(strip(a)==strip(b) for a,b in zip(r,base));runs.append({'trial':trial,'method':method,'wall_seconds':elapsed,'returncode':p.returncode,'records':len(r),'semantic_match':equal,'stderr':p.stderr});(OUT/'runs.json').write_text(json.dumps(runs,indent=2)+'\n');print('trial',trial,method,'seconds',round(elapsed,3),'match',equal,flush=True)
legacy=statistics.median(r['wall_seconds'] for r in runs if r['method']=='legacy');journal=statistics.median(r['wall_seconds'] for r in runs if r['method']=='journal');pairs=[]
for i in range(1,4):
 a=next(r for r in runs if r['trial']==i and r['method']=='legacy');b=next(r for r in runs if r['trial']==i and r['method']=='journal');pairs.append({'trial':i,'ratio':b['wall_seconds']/a['wall_seconds']})
summary={'runs':len(runs),'completed':sum(r['returncode']==0 for r in runs),'semantic_matches':sum(r['semantic_match'] for r in runs),'legacy_median_seconds':legacy,'journal_median_seconds':journal,'paired_ratios':pairs,'median_paired_ratio':statistics.median(x['ratio'] for x in pairs),'timing_scope':'complete Python subprocess including import, initialization, calculations and result writes; excludes creation of cloned workdir; same machine only, 3 repetitions','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps(summary),flush=True)
if summary['completed']!=6 or summary['semantic_matches']!=6:raise SystemExit(1)
