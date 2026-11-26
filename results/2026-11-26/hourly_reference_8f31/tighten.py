import sys; sys.path.insert(0, '/workspace/scratch/aa60c3a9e290/neural-pde-sync-20261112/scripts')
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,importlib.util
D=Path(__file__).resolve().parent;R=Path('/workspace/scratch/aa60c3a9e290/neural-pde-sync-20261112');P=Path('/workspace/scratch/887c25a726e2/round_20261126_2300_VkzGGB');c=json.loads((D/'config.json').read_text());sp=importlib.util.spec_from_file_location('ref321',R/'scripts/321_fixed_testset_references_20261112.py');m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
fixture=P/'inputs/frozen_networks.json';nets=json.loads(fixture.read_text())['networks'];n=next(x for x in nets if x['label']==c['network']);r=m.PolyReference(n).maximum(relative=F(c['relative']),absolute=F(c['absolute']),cap=c['cap']);(D/'reference.json').write_text(json.dumps(r,indent=2)+'\n');print('converged',r['converged'],'splits',r['splits'],flush=True)
