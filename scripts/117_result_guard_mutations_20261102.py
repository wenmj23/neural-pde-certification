"""Exercise script 116 on copied full-domain results; never alters original data."""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys,tempfile
GUARD=Path(__file__).with_name('116_verify_research_round_20261102.py')
def main():
 p=argparse.ArgumentParser();p.add_argument('--directory',required=True,type=Path);p.add_argument('--fixtures',nargs='+',required=True,type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 tests=[]
 with tempfile.TemporaryDirectory() as tmp:
  for case in ('original','truncated','modified_without_hash','duplicate_rehashed','wrong_label_rehashed','reversed_interval_rehashed','wrong_overlap_rehashed','wrong_summary','fabricated_broad_bounds_rehashed'):
   dest=Path(tmp)/case;shutil.copytree(a.directory,dest);rp=dest/'records.json';sp=dest/'summary.json';rows=json.loads(rp.read_text());s=json.loads(sp.read_text());rehash=False
   if case=='truncated':rows.pop()
   elif case=='modified_without_hash':rows[0]['network']='invented'
   elif case=='duplicate_rehashed':rows[-1]=rows[0];rehash=True
   elif case=='wrong_label_rehashed':rows[0]['network']='invented';rehash=True
   elif case=='reversed_interval_rehashed':rows[0]['arb_lower'],rows[0]['arb_upper']='1','-1';rehash=True
   elif case=='wrong_overlap_rehashed':rows[0]['interval_overlap']=False;rehash=True
   elif case=='wrong_summary':s['midpoint_in_arb']-=1
   elif case=='fabricated_broad_bounds_rehashed':rows[0]['arb_lower'],rows[0]['arb_upper']='-1000000000000000000','1000000000000000000';rehash=True
   if case not in ('original','wrong_summary'):
    rp.write_text(json.dumps(rows,indent=2)+'\n')
    if rehash:s['records_sha256']=hashlib.sha256(rp.read_bytes()).hexdigest()
   if case=='wrong_summary' or rehash:sp.write_text(json.dumps(s,indent=2)+'\n')
   cmd=[sys.executable,str(GUARD),'full',str(dest),'--fixtures',*[str(x) for x in a.fixtures]];r=subprocess.run(cmd,capture_output=True,text=True);expected=0 if case in ('original','fabricated_broad_bounds_rehashed') else 1
   tests.append({'case':case,'expected_exit':expected,'actual_exit':r.returncode,'as_expected':r.returncode==expected,'output':json.loads(r.stdout)})
 result={'checks':tests,'passed':sum(t['as_expected'] for t in tests),'expected':len(tests),'scope':'Fabricated yet internally consistent broad bounds are intentionally accepted; hashes and consistency are not evidence of computation authenticity.'}
 (a.out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));return int(result['passed']!=result['expected'])
if __name__=='__main__':sys.exit(main())
