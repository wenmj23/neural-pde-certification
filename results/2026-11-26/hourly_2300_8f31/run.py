from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib,importlib.util,time,itertools
R=Path('/workspace/scratch/aa60c3a9e290/neural-pde-sync-20261112');P=Path('/workspace/scratch/887c25a726e2/round_20261126_1800_CsqAVy');D=Path(__file__).parent;sys.path.insert(0,str(R/'scripts'));sp=importlib.util.spec_from_file_location('v320',R/'scripts/320_evaluate_fixed_testset_20261112.py');m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();c=json.loads((D/'protocol.json').read_text());f=P/'training/frozen_networks.json';ts=json.loads((P/'training/summary.json').read_text());assert sha(f)==ts['file_sha256']['frozen_networks.json'];nets=json.loads(f.read_text())['networks'];rows=[]
for n in nets:
 for degree,mode,budget in itertools.product(c['degrees'],c['mesh_modes'],c['budgets']):
  start=time.perf_counter();err=None
  try:
   p,v=m.construct(f,sha(f),n,degree,1,c);v.uniform(16)
   if mode=='uniform':v.uniform(budget)
   else:v.refine_to(budget)
   b=m.bounds(v.certify());nodes=[str(x) for x in v.nodes];row={'network':n['label'],'source':n['source_name'],'degree':degree,'mode':mode,'budget':budget,'status':'finite',**b,'nodes':nodes,'residual_evaluations':v.residual_evaluations,'derivative_interval_requests':v.curvature_evaluations,'underlying_derivative_calls':p.underlying_derivative_calls}
  except Exception as e:row={'network':n['label'],'source':n['source_name'],'degree':degree,'mode':mode,'budget':budget,'status':'error','error':repr(e)}
  row['seconds']=time.perf_counter()-start;rows.append(row)
 (D/('case_'+n['label']+'.json')).write_text(json.dumps([x for x in rows if x['network']==n['label']],indent=2)+'\n');print(n['label'],flush=True)
(D/'records.json').write_text(json.dumps(rows,indent=2)+'\n');sources=['320_evaluate_fixed_testset_20261112.py','arb_polynomial_forcing.py','225_derivative_split_20261105.py','151_new_pinn_matched_baselines_20261102.py','arb_cubic_reconstruction.py','arb_shaped_remainder.py','arb_deep_network.py','arb_reconstruction.py','arb_certification.py','arb_problem_families.py','interval_certification.py'];s={'expected':c['expected'],'records':len(rows),'errors':sum(x['status']=='error' for x in rows),'fixture_sha256':sha(f),'protocol_sha256':sha(D/'protocol.json'),'run_sha256':sha(D/'run.py'),'records_sha256':sha(D/'records.json'),'sources':{n:sha(R/'scripts'/n) for n in sources},'cases':{p.name:sha(p) for p in D.glob('case_*.json')}};(D/'summary.json').write_text(json.dumps(s,indent=2)+'\n');print(json.dumps({k:v for k,v in s.items() if k not in ['sources','cases']}))
