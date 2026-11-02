"""Exercise script 123 on temporary copies, including its known blind spot."""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys,tempfile
GUARD=Path(__file__).with_name('123_verify_multiple_rounds_20261102.py')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--fixtures',nargs='+',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 cases=['original','truncated','wrong_source_rehashed','reference_changed_no_hash','wrong_scale_rehashed','wrong_tight_parent_rehashed','wrong_decision_rehashed','wrong_mp_value_rehashed','fabricated_wide_arb_rehashed'];results=[]
 with tempfile.TemporaryDirectory() as temp:
  for case in cases:
   root=Path(temp)/case;shutil.copytree(a.root,root)
   folder={'truncated':'transferred_forcing','wrong_source_rehashed':'transferred_forcing','reference_changed_no_hash':'transferred_forcing','wrong_scale_rehashed':'metamorphic_scaling','wrong_tight_parent_rehashed':'tightened_reference','wrong_decision_rehashed':'decision_boundaries','wrong_mp_value_rehashed':'transferred_derivatives','fabricated_wide_arb_rehashed':'transferred_derivatives'}.get(case)
   if folder:
    d=root/folder;sp=d/'summary.json';s=json.loads(sp.read_text());rp=d/('references.json' if case=='reference_changed_no_hash' else 'records.json');rows=json.loads(rp.read_text());rehash=case.endswith('_rehashed')
    if case=='truncated':rows.pop()
    elif case=='wrong_source_rehashed':rows[0]['source']='invented'
    elif case=='reference_changed_no_hash':rows[0]['lower']='0'
    elif case=='wrong_scale_rehashed':rows[0]['normalized_upper']='0'
    elif case=='wrong_tight_parent_rehashed':rows[0]['certificate_lower']='0'
    elif case=='wrong_decision_rehashed':rows[0]['actual']='within_tolerance'
    elif case=='wrong_mp_value_rehashed':rows[0]['mp70']=str(10**100)
    elif case=='fabricated_wide_arb_rehashed':rows[0]['arb_lower'],rows[0]['arb_upper']=str(-10**100),str(10**100)
    rp.write_text(json.dumps(rows,indent=2)+'\n')
    if rehash:s['records_sha256']=hashlib.sha256(rp.read_bytes()).hexdigest();sp.write_text(json.dumps(s,indent=2)+'\n')
   cmd=[sys.executable,str(GUARD),'--root',str(root),'--fixtures',*[str(x) for x in a.fixtures]];r=subprocess.run(cmd,capture_output=True,text=True);expected=0 if case in ('original','fabricated_wide_arb_rehashed') else 1
   results.append({'case':case,'expected_exit':expected,'actual_exit':r.returncode,'as_expected':r.returncode==expected,'output':json.loads(r.stdout)})
 s={'checks':results,'passed':sum(x['as_expected'] for x in results),'expected':len(results),'scope':'Hashes and consistency do not authenticate results; deliberately fabricated broad bounds are expected to pass.'};(a.out/'summary.json').write_text(json.dumps(s,indent=2)+'\n');print(json.dumps(s));return int(s['passed']!=s['expected'])
if __name__=='__main__':sys.exit(main())
