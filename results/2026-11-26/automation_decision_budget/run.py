from pathlib import Path
from fractions import Fraction as F
import subprocess,json,itertools,hashlib,time
root=Path('/workspace/scratch/aa60c3a9e290/neural-pde-sync-20261112'); out=Path(__file__).parent
py='/workspace/scratch/aa60c3a9e290/pde-venv/bin/python'
rows=[]
for i,(name,deg,tol,budget) in enumerate(itertools.product(['shallow','deep'],[2,3],['1e-3','1e-4'],[16,32])):
 dest=out/f'case_{i:02d}'
 argv=[py,str(root/'scripts/certify.py'),'--case',str(root/f'configs/examples/{name}_sine.json'),'--out',str(dest),'--goal','decision','--degree',str(deg),'--error-tolerance',tol,'--max-intervals',str(budget)]
 p=subprocess.run(argv,capture_output=True,text=True,timeout=60,cwd=root); r=json.loads((dest/'result.json').read_text());e=r['enclosure'];lo=F(int(e['lower']['numerator']),int(e['lower']['denominator']));hi=F(int(e['upper']['numerator']),int(e['upper']['denominator'])); t=F(tol);status='within_tolerance' if hi<=t else 'exceeds_tolerance' if lo>t else 'undetermined';assert status==r['accuracy_decision'];assert p.returncode==(0 if r['goal_met'] else 2)
 rows.append({'case':name,'degree':deg,'tolerance':tol,'budget':budget,'lower':str(lo),'upper':str(hi),'decision':status,'stop':r['stop_reason'],'exit_code':p.returncode,'result_sha256':hashlib.sha256((dest/'result.json').read_bytes()).hexdigest()})
(out/'records.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps({'count':len(rows),'status':{k:sum(x['decision']==k for x in rows) for k in ['within_tolerance','exceeds_tolerance','undetermined']},'records_sha256':hashlib.sha256((out/'records.json').read_bytes()).hexdigest()}))
