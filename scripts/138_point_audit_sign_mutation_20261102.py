"""Check 134 exits nonzero on an isolated sign mutation of the source term.

Uses one supplied frozen network; the production implementation is untouched.
The output intentionally contains the negative control's audit failures.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,importlib.util,json,sys

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--fixture',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 nets=json.loads(a.fixture.read_bytes())['networks']
 if len(nets)!=1:raise ValueError('negative control requires exactly one network fixture')
 file=Path(__file__).with_name('134_point_audit_polynomial_20261102.py');spec=importlib.util.spec_from_file_location('audit134',file);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);base=m.prototype()
 class WrongSign(base):
  def source_derivative(self,x,j):return -super().source_derivative(x,j)
 m.prototype=lambda:WrongSign
 rc=m.run([a.fixture],a.out/'fault_audit','all');s=json.loads((a.out/'fault_audit/summary.json').read_text());path=a.out/'fault_audit/records.json';rows=json.loads(path.read_text())
 if sha(path)!=s['records_sha256'] or len(rows)!=48:raise ValueError('negative control record integrity')
 mismatches=[r for r in rows if not r['contained']];wrong=[r for r in rows if (not r['contained'])!=(F(r['forcing_derivative'])!=0) or not r['stable']]
 summary={'audit_return_code':rc,'expected_return_code':1,'records':len(rows),'observed_failures':len(mismatches),'expected_failures':sum(F(r['forcing_derivative'])!=0 for r in rows),'unexpected_results':len(wrong),'audit_records_sha256':s['records_sha256'],'fixture_sha256':sha(a.fixture),'test_passed':rc==1 and not wrong,'scope':'one-network intentional wrong-sign control; numerical point diagnostic only'}
 (a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary));return int(not summary['test_passed'])
if __name__=='__main__':sys.exit(main())
