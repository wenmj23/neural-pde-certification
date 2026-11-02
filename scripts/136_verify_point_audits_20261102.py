"""Post-exit manifest and exact arithmetic check for scripts 134 and 135.

Checks saved point intervals, numerical approximations and fault expectations;
does not rerun either backend or prove the reference numbers.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import argparse,hashlib,json,sys
TRANSFER={'quadratic_12x2':[0,0,12],'cubic_20x3':[0,0,0,20],'quintic_42x5':[0,0,0,0,0,42]};P=('1/7','1/2','6/7');O=(0,2,3,4)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(b,msg):
 if not b:raise ValueError(msg)
def key(r):return (r['fixture_sha256'],r['network'],r['source'],r['x'],r['order'])
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--normal',required=True,type=Path);ap.add_argument('--fault',required=True,type=Path);ap.add_argument('--fixtures',required=True,nargs='+',type=Path);a=ap.parse_args()
 fixture_nets={(sha(path),n['label']):n for path in a.fixtures for n in json.loads(path.read_bytes())['networks']};need(len(fixture_nets)>0,'no fixtures')
 normal_summary=json.loads((a.normal/'summary.json').read_text());mode=normal_summary['mode'];need(mode in ('all','original','transfer'),'unsupported mode');sources=('original',*TRANSFER) if mode=='all' else (('original',) if mode=='original' else tuple(TRANSFER));manifest={(h,n,src,x,o) for (h,n) in fixture_nets for src in sources for x in P for o in O};all_records=[]
 for directory,fault in ((a.normal,False),(a.fault,True)):
  s=json.loads((directory/'summary.json').read_text());path=directory/'records.json';need(sha(path)==s['records_sha256'],'records hash');need(sha(Path(__file__).with_name('134_point_audit_polynomial_20261102.py'))==s['script_sha256'],'audit script hash');need(sha(Path(__file__).with_name('107_polynomial_derivative_crosscheck_20261101.py'))==s['point_reference_sha256'],'reference script hash');need(s['mode']==mode and s['test_fault']==('omit_source' if fault else 'none') and s['transfer_sources']==TRANSFER and s['points']==list(P) and s['orders']==list(O),'configuration mismatch')
  records=json.loads(path.read_text());keys=[key(r) for r in records];need(len(keys)==len(manifest)==len(set(keys))==s['records']==s['expected'] and set(keys)==manifest,'manifest mismatch')
  for r in records:
   net=fixture_nets[r['fixture_sha256'],r['network']];coeff=net['source_coefficients'] if r['source']=='original' else TRANSFER[r['source']];need(r['source_coefficients']==coeff,'source coefficients mismatch');x=F(r['x']);j=r['order'];forcing=sum((F(c*factorial(i),factorial(i-j))*x**(i-j) for i,c in enumerate(coeff) if i>=j),F(0));lo,hi,low,high=map(F,(r['arb_lower'],r['arb_upper'],r['mp60'],r['mp90']));inside=lo<=high<=hi;stable=abs(low-high)<=F(1,10**50)*max(F(1),abs(high));expected=not (fault and forcing!=0)
   need(lo<=hi and F(r['forcing_derivative'])==forcing,'invalid endpoint/forcing');need(r['contained'] is inside and r['stable'] is stable and r['expected_contained'] is expected and r['passed'] is (stable and inside==expected),'stored check mismatch')
  for field,fn in [('stable',lambda r:r['stable']),('contained',lambda r:r['contained']),('expected_mismatch',lambda r:not r['expected_contained']),('passed',lambda r:r['passed'])]:need(s[field]==sum(fn(r) for r in records),'summary '+field)
  need(s['passed']==len(records) and not s['failures'],'unexpected failures');all_records.append({key(r):r for r in records})
 normal,bug=all_records
 for k in manifest:need(normal[k]['mp60']==bug[k]['mp60'] and normal[k]['mp90']==bug[k]['mp90'],'reference changed between runs')
 print(json.dumps({'ok':True,'normal_records':len(normal),'fault_records':len(bug),'fault_expected_mismatches':sum(not x['expected_contained'] for x in bug.values()),'scope':'readback validation of saved arithmetic and request coverage, not proof of numerical backends'}));return 0
if __name__=='__main__':
 try:sys.exit(main())
 except (OSError,ValueError,TypeError,KeyError,ZeroDivisionError) as e:print(json.dumps({'ok':False,'error':str(e)}));sys.exit(1)
