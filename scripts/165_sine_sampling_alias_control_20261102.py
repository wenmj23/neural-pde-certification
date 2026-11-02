"""Demonstrate finite-point blindness with intentionally chosen k=7 and k=77.

Both networks reuse one frozen k=1 network's weights with only frequency and
label changed. They are artificial controls, not trained k=7/k=77 PINNs.
"""
from copy import deepcopy
from pathlib import Path
import argparse,hashlib,importlib.util,json,os,sys

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2,allow_nan=False)+'\n').encode();tmp=p.with_suffix('.tmp')
 with tmp.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--training',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False);source=a.training/'frozen_networks.json';s=json.loads((a.training/'summary.json').read_text())
 if sha(source)!=s['fixture_sha256']:raise ValueError('training fixture hash')
 nets=json.loads(source.read_bytes())['networks'];matches=[n for n in nets if n['label']=='sine_k1_depth1_seed501']
 if len(matches)!=1:raise ValueError('expected source model missing')
 derived=[]
 for k in (7,77):
  n=deepcopy(matches[0]);n['label']=f'control_untrained_k{k}_weights_from_k1';n['frequency']=k;derived.append(n)
 fixture=a.out/'control_fixture.json';atom(fixture,{'origin_sha256':sha(source),'note':'artificial frequency-mismatched controls, not trained k7/k77','networks':derived})
 file=Path(__file__).with_name('163_sine_derivative_point_audit_20261102.py');spec=importlib.util.spec_from_file_location('audit163',file);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 normal_rc=m.run([fixture],a.out/'normal',False);fault_rc=m.run([fixture],a.out/'omission',True);normal=json.loads((a.out/'normal/summary.json').read_text());fault=json.loads((a.out/'omission/summary.json').read_text());rows=json.loads((a.out/'omission/records.json').read_text());counts={}
 for k in (7,77):
  for group,points in (('three_original',('1/7','2/7','5/7')),('extra_point',('1/11',))):
   subset=[r for r in rows if r['frequency']==k and r['x'] in points];counts[f'{k}/{group}']={'total':len(subset),'undetected':sum(r['contained'] for r in subset),'detected':sum(not r['contained'] for r in subset),'source_zero_exact':sum(r['source_zero_exact'] for r in subset)}
 expected={'7/three_original':{'total':12,'undetected':9,'detected':3,'source_zero_exact':9},'7/extra_point':{'total':4,'undetected':0,'detected':4,'source_zero_exact':0},'77/three_original':{'total':12,'undetected':9,'detected':3,'source_zero_exact':9},'77/extra_point':{'total':4,'undetected':3,'detected':1,'source_zero_exact':3}}
 result={'original_fixture_sha256':sha(source),'derived_fixture_sha256':sha(fixture),'normal_records_sha256':normal['records_sha256'],'omission_records_sha256':fault['records_sha256'],'script_sha256':sha(Path(__file__)),'counts':counts,'expected_counts':expected,'normal_rc':normal_rc,'omission_rc':fault_rc,'test_passed':normal_rc==0 and fault_rc==0 and counts==expected,'scope':'constructed aliasing counterexample; intentionally faulty subclass and synthetic frequency, not a solver failure'};atom(a.out/'summary.json',result);print(json.dumps(result),flush=True);return int(not result['test_passed'])
if __name__=='__main__':sys.exit(main())
