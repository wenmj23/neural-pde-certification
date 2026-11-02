"""Record and compare the initial 60/90 and rerun 80/110 precision gates."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,os,sys

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2)+'\n').encode();tmp=p.with_suffix('.tmp')
 with tmp.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def key(r):return (r['fixture_sha256'],r['network'],r['x'],r['order'])
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--initial',required=True,type=Path);ap.add_argument('--recheck',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 old_s=json.loads((a.initial/'summary.json').read_text());new_s=json.loads((a.recheck/'summary.json').read_text());op=a.initial/'records.json';np=a.recheck/'records.json'
 if sha(op)!=old_s['records_sha256'] or sha(np)!=new_s['records_sha256']:raise ValueError('record hash')
 old=json.loads(op.read_text());new=json.loads(np.read_text());o={key(r):r for r in old};n={key(r):r for r in new}
 if len(o)!=len(old) or len(n)!=len(new) or set(o)!=set(n):raise ValueError('manifest mismatch')
 rows=[]
 for k,r in o.items():
  q=n[k];l,h=F(r['mp60']),F(r['mp90']);l2,h2=F(q['mp80']),F(q['mp110']);old_rel=abs(l-h)/max(F(1),abs(h));new_rel=abs(l2-h2)/max(F(1),abs(h2));high_rel=abs(h-h2)/max(F(1),abs(h2));old_stable=old_rel<=F(1,10**50);new_stable=new_rel<=F(1,10**50)
  if r['stable']!=old_stable or q['stable']!=new_stable:raise ValueError('stability mismatch')
  if not old_stable:rows.append({'network':r['network'],'frequency':r['frequency'],'x':r['x'],'order':r['order'],'old_relative_difference':str(old_rel),'new_relative_difference':str(new_rel),'old_high_vs_new_high_relative_difference':str(high_rel),'old_stable':old_stable,'new_stable':new_stable,'both_point_contained':bool(r['contained'] and q['contained'])})
 atom(a.out/'records.json',rows);s={'initial_records_sha256':old_s['records_sha256'],'recheck_records_sha256':new_s['records_sha256'],'initial_audit_script_sha256':old_s['script_sha256'],'recheck_audit_script_sha256':new_s['script_sha256'],'initial_digits':[60,90],'recheck_digits':[80,110],'threshold':'1e-50 * max(1, abs(high-precision value))','initial_records':len(old),'initial_unstable':len(rows),'recovered':sum(r['new_stable'] for r in rows),'both_contained':sum(r['both_point_contained'] for r in rows),'records_sha256':sha(a.out/'records.json'),'script_sha256':sha(Path(__file__)),'scope':'numerical stability incident, not proof or evidence of a certificate failure'};atom(a.out/'summary.json',s);print(json.dumps(s),flush=True);return int(len(rows)!=4 or s['recovered']!=4)
if __name__=='__main__':sys.exit(main())
