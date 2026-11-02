"""Post-hoc near-threshold classification stress on script 139 outputs.

Thresholds are derived from prior reference midpoint and thus are not blind.
No new certificates are calculated and references never select a certifier mesh.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json,os,sys
DELTAS=(F(1,1000),F(1,10000),F(1,100000));SIGNS=(-1,1)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2)+'\n').encode();tmp=p.with_suffix('.tmp')
 with tmp.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def classify(l,u,t):return 'within' if u<=t else ('exceeds' if l>t else 'undetermined')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False);s=json.loads((a.input/'summary.json').read_text());p=a.input/'records.json'
 if sha(p)!=s['records_sha256']:raise ValueError('source hash mismatch')
 inputs=json.loads(p.read_text());rows=[]
 for r in inputs:
  lo,hi,rl,ru=map(F,(r['lower'],r['upper'],r['reference_lower'],r['reference_upper']));mid=(rl+ru)/2
  for delta in DELTAS:
   for sign in SIGNS:
    t=mid*(1+sign*delta);d,rd=classify(lo,hi,t),classify(rl,ru,t);compatible=d=='undetermined' or rd=='undetermined' or d==rd
    rows.append({'fixture_sha256':r['fixture_sha256'],'network':r['network'],'source':r['source'],'degree':r['degree'],'intervals':r['intervals'],'relative_offset':str(delta),'side':sign,'tolerance':str(t),'decision':d,'reference_decision':rd,'compatible':compatible})
 atom(a.out/'records.json',rows);summary={'input_sha256':s['records_sha256'],'input_records':len(inputs),'relative_offsets':[str(x) for x in DELTAS],'sides':SIGNS,'expected':len(inputs)*len(DELTAS)*len(SIGNS),'records':len(rows),'records_sha256':sha(a.out/'records.json'),'script_sha256':sha(Path(__file__)),'counts':{},'reference_undetermined':sum(r['reference_decision']=='undetermined' for r in rows),'incompatible':sum(not r['compatible'] for r in rows),'scope':'post-hoc thresholds from reference midpoint; diagnostic stress only, not blinded'}
 for d in (2,3):
  for n in (8,16,32,64):
   for delta in DELTAS:
    group=[r for r in rows if r['degree']==d and r['intervals']==n and r['relative_offset']==str(delta)];summary['counts'][f'{d}/{n}/{delta}']={'within':sum(r['decision']=='within' for r in group),'exceeds':sum(r['decision']=='exceeds' for r in group),'undetermined':sum(r['decision']=='undetermined' for r in group),'reference_undetermined':sum(r['reference_decision']=='undetermined' for r in group),'incompatible':sum(not r['compatible'] for r in group)}
 if len(rows)!=summary['expected']:raise ValueError('count mismatch')
 atom(a.out/'summary.json',summary);print(json.dumps(summary),flush=True);return int(summary['incompatible']!=0)
if __name__=='__main__':sys.exit(main())
