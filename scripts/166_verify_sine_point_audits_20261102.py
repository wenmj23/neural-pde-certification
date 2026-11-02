"""Independently read back sine point audits and deliberate alias control."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from copy import deepcopy
import argparse,hashlib,json,sys
P=('1/7','2/7','5/7','1/11');O=(0,2,3,4)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(x,msg):
 if not x:raise ValueError(msg)
def key(r):return (r['fixture_sha256'],r['network'],r['x'],r['order'])
def source_zero(k,x,j):return (k*x.numerator)%x.denominator==0 if j%2==0 else (2*k*x.numerator)%(2*x.denominator)==x.denominator
def verify(normal,fault,fixture):
 nets=json.loads(fixture.read_bytes())['networks'];h=sha(fixture);bylabel={n['label']:n for n in nets};need(len(bylabel)==len(nets)>0,'networks');expected=set(product((h,),bylabel,P,O));results=[]
 for directory,isfault in ((normal,False),(fault,True)):
  s=json.loads((directory/'summary.json').read_text());p=directory/'records.json';need(sha(p)==s['records_sha256'] and sha(Path(__file__).with_name('163_sine_derivative_point_audit_20261102.py'))==s['script_sha256'],'record/script hash');need(s['inputs']==[{'sha256':h,'networks':len(nets)}] and s['points']==list(P) and s['orders']==list(O) and s['reference_digits']==[80,110] and s['test_fault']==('omit_forcing' if isfault else 'none'),'config')
  rows=json.loads(p.read_text());keys=[key(r) for r in rows];need(len(rows)==len(set(keys))==len(expected)==s['records']==s['expected'] and set(keys)==expected,'manifest')
  for r in rows:
   k=bylabel[r['network']]['frequency'];x=F(r['x']);order=r['order'];need(r['frequency']==k,'frequency');zero=source_zero(k,x,order);lo,hi,low,high=map(F,(r['arb_lower'],r['arb_upper'],r['mp80'],r['mp110']));stable=abs(low-high)<=F(1,10**50)*max(F(1),abs(high));inside=lo<=high<=hi;expected_inside=not isfault or zero
   need(lo<=hi and r['source_zero_exact'] is zero and r['stable'] is stable and r['contained'] is inside and r['expected_contained'] is expected_inside and r['passed'] is (stable and inside==expected_inside),'record arithmetic')
  for field,fn in [('stable',lambda r:r['stable']),('contained',lambda r:r['contained']),('source_zero_exact',lambda r:r['source_zero_exact']),('expected_mismatch',lambda r:not r['expected_contained']),('passed',lambda r:r['passed'])]:need(s[field]==sum(fn(r) for r in rows),'summary '+field)
  need(s['passed']==len(rows) and not s['failures'],'audit reported failure');results.append({key(r):r for r in rows})
 normalrows,faultrows=results
 for k in expected:need(normalrows[k]['mp80']==faultrows[k]['mp80'] and normalrows[k]['mp110']==faultrows[k]['mp110'],'references differ across normal and fault')
 return normalrows,faultrows
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--training',type=Path,required=True);ap.add_argument('--normal',type=Path,required=True);ap.add_argument('--fault',type=Path,required=True);ap.add_argument('--alias',type=Path,required=True);ap.add_argument('--alias-fixture',type=Path);a=ap.parse_args();fixture=a.training/'frozen_networks.json';ts=json.loads((a.training/'summary.json').read_text());need(sha(fixture)==ts['fixture_sha256'],'training fixture');normal,fault=verify(a.normal,a.fault,fixture)
 p=a.alias;alias=json.loads((p/'summary.json').read_text());derived=a.alias_fixture or (p/'control_fixture.json');need(sha(fixture)==alias['original_fixture_sha256'] and sha(derived)==alias['derived_fixture_sha256'] and sha(Path(__file__).with_name('165_sine_sampling_alias_control_20261102.py'))==alias['script_sha256'],'alias fixture/script');original=[n for n in json.loads(fixture.read_bytes())['networks'] if n['label']=='sine_k1_depth1_seed501'];need(len(original)==1,'original network');nets=json.loads(derived.read_bytes())['networks'];need(len(nets)==2,'derived networks')
 for k in (7,77):
  n=deepcopy(original[0]);n['label']=f'control_untrained_k{k}_weights_from_k1';n['frequency']=k;need(sum(x==n for x in nets)==1,'derived weights')
 an,af=verify(p/'normal',p/'omission',derived);need(alias['normal_records_sha256']==sha(p/'normal/records.json') and alias['omission_records_sha256']==sha(p/'omission/records.json'),'alias records');counts={}
 for k in (7,77):
  for group,points in (('three_original',('1/7','2/7','5/7')),('extra_point',('1/11',))):
   subset=[r for r in af.values() if r['frequency']==k and r['x'] in points];counts[f'{k}/{group}']={'total':len(subset),'undetected':sum(r['contained'] for r in subset),'detected':sum(not r['contained'] for r in subset),'source_zero_exact':sum(r['source_zero_exact'] for r in subset)}
 expected_counts={'7/three_original':{'total':12,'undetected':9,'detected':3,'source_zero_exact':9},'7/extra_point':{'total':4,'undetected':0,'detected':4,'source_zero_exact':0},'77/three_original':{'total':12,'undetected':9,'detected':3,'source_zero_exact':9},'77/extra_point':{'total':4,'undetected':3,'detected':1,'source_zero_exact':3}}
 need(alias['test_passed'] is True and alias['normal_rc']==alias['omission_rc']==0 and alias['counts']==counts==expected_counts and alias['expected_counts']==expected_counts,'alias counts');print(json.dumps({'ok':True,'normal_records':len(normal),'fault_records':len(fault),'alias_normal':len(an),'alias_fault':len(af),'alias_counts':counts,'scope':'readback consistency and exact zero arithmetic; not proof of either numerical backend'}));return 0
if __name__=='__main__':
 try:sys.exit(main())
 except (OSError,ValueError,KeyError,TypeError,ZeroDivisionError) as e:print(json.dumps({'ok':False,'error':str(e)}));sys.exit(1)
