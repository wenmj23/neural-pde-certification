"""Run and stress-test resumable decision persistence on fixed frozen inputs."""
from pathlib import Path
import hashlib,importlib.util,json,shutil,subprocess,sys
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261101/resumable_audit';OUT.mkdir(parents=True,exist_ok=False)
runner=ROOT/'scripts/85_resumable_decisions_20261101.py';checks=[];events=[]
def call(name,args,expected):
 p=subprocess.run([sys.executable,str(runner),*args],cwd=ROOT,capture_output=True,text=True,timeout=180);events.append({'name':name,'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr});check(name+'_exit',p.returncode==expected);return p
def check(name,ok):
 checks.append({'name':name,'passed':bool(ok)});(OUT/'checks.json').write_text(json.dumps(checks,indent=2)+'\n');print(name,ok,flush=True)
normal=OUT/'normal';partial=OUT/'interrupted';bad_tail=OUT/'bad_tail';dup=OUT/'duplicate';tamper=OUT/'tampered_export'
call('normal_run',['--out',str(normal)],0)
call('controlled_stop',['--out',str(partial),'--simulate-stop-after','17'],3)
before=(partial/'records.jsonl').read_bytes();check('stopped_with_17_rows',len(before.splitlines())==17);check('no_final_summary_during_stop',not (partial/'summary.json').exists())
shutil.copytree(partial,bad_tail);t=bad_tail/'records.jsonl';t.write_bytes(t.read_bytes()[:-1]);digest=hashlib.sha256(t.read_bytes()).hexdigest();z=call('reject_incomplete_line',['--out',str(bad_tail),'--resume'],1);check('incomplete_line_left_unchanged',hashlib.sha256(t.read_bytes()).hexdigest()==digest);check('incomplete_line_reason','incomplete last line' in z.stderr)
shutil.copytree(partial,dup);q=dup/'records.jsonl';lines=q.read_bytes().splitlines(keepends=True);q.write_bytes(q.read_bytes()+lines[0]);digest=hashlib.sha256(q.read_bytes()).hexdigest();z=call('reject_duplicate',['--out',str(dup),'--resume'],1);check('duplicate_left_unchanged',hashlib.sha256(q.read_bytes()).hexdigest()==digest);check('duplicate_reason','Duplicate request' in z.stderr)
call('resume_after_stop',['--out',str(partial),'--resume'],0);after=(partial/'records.jsonl').read_bytes();check('resumed_kept_first_17_rows',after.startswith(before));check('resumed_has_72_rows',len(after.splitlines())==72)
spec=importlib.util.spec_from_file_location('journal_runner',runner);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
r1,s1=m.verify_export(normal);r2,s2=m.verify_export(partial);base=json.loads((ROOT.parent/'run_20261101/poly_fixed_tolerance/records.json').read_text());strip=lambda x:{k:v for k,v in x.items() if k!='elapsed_seconds'}
check('normal_matches_existing_all_72',len(r1)==len(base)==72 and all(strip(a)==strip(b) for a,b in zip(r1,base)));check('resumed_matches_existing_all_72',len(r2)==len(base)==72 and all(strip(a)==strip(b) for a,b in zip(r2,base)));check('no_reference_disagreements',all(r['compatible'] for r in r1+r2))
shutil.copytree(normal,tamper);f=tamper/'records.json';f.write_bytes(f.read_bytes()+b' ');bad=False
try:m.verify_export(tamper)
except ValueError as e:bad='checksum' in str(e)
check('modified_export_detected',bad)
for p in (normal,partial):
 rows=json.loads((p/'records.json').read_text());summary=json.loads((p/'summary.json').read_text());check(p.name+'_summary_count',len(rows)==summary['expected_count']==summary['records_count']==72)
(OUT/'events.json').write_text(json.dumps(events,indent=2)+'\n');summary={'checks':len(checks),'passed':sum(x['passed'] for x in checks),'failed':[x for x in checks if not x['passed']],'runner_sha256':hashlib.sha256(runner.read_bytes()).hexdigest(),'test_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'note':'Controlled stop follows a fully flushed record; no proof for every abrupt hardware or power failure; original anomaly remains unexplained.'};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps(summary),flush=True)
if summary['failed']:raise SystemExit(1)
