"""Exercise the point-audit result guard on temporary copied data."""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys,tempfile
GUARD=Path(__file__).with_name('136_verify_point_audits_20261102.py')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--normal',type=Path,required=True);ap.add_argument('--fault',type=Path,required=True);ap.add_argument('--fixtures',type=Path,nargs='+',required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 cases=['original','truncated','duplicate_rehashed','wrong_source_rehashed','wrong_reference_rehashed','wrong_flag_rehashed','wrong_summary','fabricated_wide_bounds_rehashed'];checks=[]
 with tempfile.TemporaryDirectory() as tmp:
  for i,case in enumerate(cases):
   root=Path(tmp)/str(i);root.mkdir();normal=root/'normal';fault=root/'fault';shutil.copytree(a.normal,normal);shutil.copytree(a.fault,fault);rp=normal/'records.json';sp=normal/'summary.json';rows=json.loads(rp.read_text());summary=json.loads(sp.read_text());rehash=case.endswith('_rehashed')
   if case=='truncated':rows.pop()
   elif case=='duplicate_rehashed':rows[-1]=rows[0]
   elif case=='wrong_source_rehashed':rows[0]['source']='unknown'
   elif case=='wrong_reference_rehashed':rows[0]['mp90']=str(10**100)
   elif case=='wrong_flag_rehashed':rows[0]['contained']=False
   elif case=='wrong_summary':summary['passed']-=1
   elif case=='fabricated_wide_bounds_rehashed':rows[0]['arb_lower'],rows[0]['arb_upper']=str(-10**100),str(10**100)
   if case not in ('original','wrong_summary'):
    rp.write_text(json.dumps(rows,indent=2)+'\n')
    if rehash:summary['records_sha256']=hashlib.sha256(rp.read_bytes()).hexdigest()
   if rehash or case=='wrong_summary':sp.write_text(json.dumps(summary,indent=2)+'\n')
   cmd=[sys.executable,str(GUARD),'--normal',str(normal),'--fault',str(fault),'--fixtures',*[str(f) for f in a.fixtures]];r=subprocess.run(cmd,capture_output=True,text=True);expected=0 if case in ('original','fabricated_wide_bounds_rehashed') else 1
   checks.append({'case':case,'expected_exit':expected,'actual_exit':r.returncode,'as_expected':r.returncode==expected,'output':json.loads(r.stdout)})
 summary={'expected':len(checks),'passed':sum(r['as_expected'] for r in checks),'checks':checks,'scope':'tests result guard behavior; fabricated consistent intervals and recomputed hashes are intentionally not detected'};(a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary));return int(summary['passed']!=summary['expected'])
if __name__=='__main__':sys.exit(main())
