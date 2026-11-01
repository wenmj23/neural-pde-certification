"""Mutation tests for script 108 on temporary copies, including a known limitation."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys,tempfile,shutil
P=Path(__file__).resolve().parent/'108_verify_result_bundle_20261101.py'
ap=argparse.ArgumentParser();ap.add_argument('--certificate',required=True,type=Path);ap.add_argument('--derivative',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);ap.add_argument('--fixture',required=True,type=Path);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
def call(p,fixture=None):
 cmd=[sys.executable,str(P),str(p)]
 if fixture:cmd.extend(['--fixtures',str(fixture)])
 r=subprocess.run(cmd,capture_output=True,text=True);return r.returncode,json.loads(r.stdout)
def recalc(p,rows):
 b=(json.dumps(rows,indent=2,ensure_ascii=False)+'\n').encode();(p/'records.json').write_bytes(b);s=json.loads((p/'summary.json').read_text());s['records_sha256']=hashlib.sha256(b).hexdigest();(p/'summary.json').write_text(json.dumps(s))
checks=[]
with tempfile.TemporaryDirectory() as tmp:
 root=Path(tmp)
 def trial(name,src,mutation,want):
  p=root/name;shutil.copytree(src,p);mutation(p);code,output=call(p);checks.append({'case':name,'expected_exit':want,'actual_exit':code,'detected_as_expected':code==want,'errors':output['errors']})
 def remove_some(p):
  f=p/'records.json';x=json.loads(f.read_text());f.write_text(json.dumps(x[:-12]))
 def duplicate_rehash(p):
  f=p/'records.json';x=json.loads(f.read_text());x[-1]=x[0];recalc(p,x)
 def wrong_count(p):
  f=p/'summary.json';s=json.loads(f.read_text());s['records']-=1;f.write_text(json.dumps(s))
 def corrupt_json(p):(p/'records.json').write_text('{oops')
 def missing(p):(p/'summary.json').unlink()
 def change_expected(p):
  f=p/'summary.json';s=json.loads(f.read_text());s['expected']-=1;f.write_text(json.dumps(s))
 def wrong_interval(p):
  f=p/'records.json';x=json.loads(f.read_text());x[0]['intervals']+=1;recalc(p,x)
 def false_total(p):
  f=p/'summary.json';s=json.loads(f.read_text());s['contained']-=1;f.write_text(json.dumps(s))
 def replace_one_label(p):
  f=p/'records.json';x=json.loads(f.read_text());old=x[0]['network']
  for z in x:
   if z['network']==old:z['network']='made_up_network'
  recalc(p,x)
 for name,src,mut,expected in [('certificate_valid',a.certificate,lambda p:None,0),('derivative_valid',a.derivative,lambda p:None,0),('truncated_derivative',a.derivative,remove_some,1),('duplicate_certificate_rehashed',a.certificate,duplicate_rehash,1),('duplicate_derivative_rehashed',a.derivative,duplicate_rehash,1),('wrong_count',a.certificate,wrong_count,1),('corrupt_json',a.certificate,corrupt_json,1),('missing_summary',a.certificate,missing,1),('wrong_expected',a.certificate,change_expected,1),('wrong_interval_rehashed',a.certificate,wrong_interval,1),('wrong_containment_total',a.derivative,false_total,1),('known_gap_renamed_network_rehashed',a.certificate,replace_one_label,0)]:trial(name,src,mut,expected)
with tempfile.TemporaryDirectory() as tmp:
 p=Path(tmp)/'spoof';shutil.copytree(a.certificate,p);x=json.loads((p/'records.json').read_text());old=x[0]['network']
 for z in x:
  if z['network']==old:z['network']='made_up_network'
 recalc(p,x);code,output=call(p,a.fixture);checks.append({'case':'renamed_network_with_trusted_fixture','expected_exit':1,'actual_exit':code,'detected_as_expected':code==1,'errors':output['errors']})
 code,output=call(a.certificate,a.fixture);checks.append({'case':'valid_with_trusted_fixture','expected_exit':0,'actual_exit':code,'detected_as_expected':code==0,'errors':output['errors']})
result={'tests':len(checks),'passed':sum(x['detected_as_expected'] for x in checks),'checks':checks,'limitation':'With no trusted fixture, self-consistent replacement of a whole network is undetectable; optional fixture verifies labels and its hash. SHA-256 alone is not authentication, and this is not a mathematical proof.'};(a.out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));sys.exit(0 if result['passed']==result['tests'] else 1)
