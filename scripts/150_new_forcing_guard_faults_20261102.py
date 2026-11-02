"""Test the readback guard against deliberate temporary-file mutations."""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys,tempfile
GUARD=Path(__file__).with_name('145_verify_new_forcing_20261102.py')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--training',required=True,type=Path);ap.add_argument('--certificates',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 checks=[];cases=('original','truncated','duplicate_rehashed','wrong_decision_rehashed','reference_changed','weights_changed')
 with tempfile.TemporaryDirectory() as td:
  for i,case in enumerate(cases):
   root=Path(td)/str(i);root.mkdir();training=root/'training';cert=root/'cert';shutil.copytree(a.training,training);shutil.copytree(a.certificates,cert)
   if case in ('truncated','duplicate_rehashed','wrong_decision_rehashed'):
    rp=cert/'records.json';sp=cert/'summary.json';rows=json.loads(rp.read_text());s=json.loads(sp.read_text())
    if case=='truncated':rows.pop()
    elif case=='duplicate_rehashed':rows[-1]=rows[0]
    else:rows[0]['decisions']['1/1000']='exceeds'
    rp.write_text(json.dumps(rows,indent=2)+'\n')
    if case!='truncated':s['records_sha256']=hashlib.sha256(rp.read_bytes()).hexdigest();sp.write_text(json.dumps(s,indent=2)+'\n')
   elif case=='reference_changed':
    p=cert/'references.json';x=json.loads(p.read_text());x[0]['lower']='0';p.write_text(json.dumps(x,indent=2)+'\n')
   elif case=='weights_changed':
    p=training/'frozen_networks.json';x=json.loads(p.read_text());x['networks'][0]['layers'][-1]['bias'][0]+=1;p.write_text(json.dumps(x,indent=2)+'\n')
   r=subprocess.run([sys.executable,str(GUARD),'--training',str(training),'--certificates',str(cert)],capture_output=True,text=True);expect=0 if case=='original' else 1;checks.append({'case':case,'expected':expect,'actual':r.returncode,'as_expected':r.returncode==expect,'output':json.loads(r.stdout)})
 result={'tests':len(checks),'passed':sum(c['as_expected'] for c in checks),'checks':checks,'scope':'controlled checks of selected accidental changes; self-consistent fabricated data and hash replacement are not generally detectable'};(a.out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));return int(result['tests']!=result['passed'])
if __name__=='__main__':sys.exit(main())
