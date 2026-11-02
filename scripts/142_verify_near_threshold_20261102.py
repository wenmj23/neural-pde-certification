"""Readback and exact decision check for the post-hoc threshold stress."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json,sys
D=(F(1,1000),F(1,10000),F(1,100000));SIDES=(-1,1)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def key(r):return (r['fixture_sha256'],r['network'],r['source'],r['degree'],r['intervals'])
def decide(l,u,t):return 'within' if u<=t else ('exceeds' if l>t else 'undetermined')
def need(v,msg):
 if not v:raise ValueError(msg)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True,type=Path);ap.add_argument('--directory',required=True,type=Path);a=ap.parse_args();s=json.loads((a.directory/'summary.json').read_text());prev=json.loads((a.input/'summary.json').read_text());need(sha(a.input/'records.json')==prev['records_sha256']==s['input_sha256'],'input hash');need(sha(a.directory/'records.json')==s['records_sha256'],'records hash');need(sha(Path(__file__).with_name('141_near_threshold_stress_20261102.py'))==s['script_sha256'],'script hash');need(s['relative_offsets']==[str(d) for d in D] and s['sides']==list(SIDES),'settings')
 parent=json.loads((a.input/'records.json').read_text());pmap={key(r):r for r in parent};need(len(pmap)==len(parent)==s['input_records'],'input keys');rows=json.loads((a.directory/'records.json').read_text());keys=[(*key(r),r['relative_offset'],r['side']) for r in rows];expected={(*k,str(d),sign) for k in pmap for d in D for sign in SIDES};need(len(rows)==len(set(keys))==len(expected)==s['records']==s['expected'] and set(keys)==expected,'record coverage')
 for r in rows:
  orig=pmap[key(r)];lo,hi,rl,ru=map(F,(orig['lower'],orig['upper'],orig['reference_lower'],orig['reference_upper']));delta=F(r['relative_offset']);t=(rl+ru)/2*(1+r['side']*delta);actual,ref=decide(lo,hi,t),decide(rl,ru,t);compatible=actual=='undetermined' or ref=='undetermined' or actual==ref;need(F(r['tolerance'])==t and r['decision']==actual and r['reference_decision']==ref and r['compatible']==compatible and compatible,'decision mismatch')
 for degree in (2,3):
  for n in (8,16,32,64):
   for d in D:
    group=[r for r in rows if r['degree']==degree and r['intervals']==n and r['relative_offset']==str(d)];counts={'within':sum(r['decision']=='within' for r in group),'exceeds':sum(r['decision']=='exceeds' for r in group),'undetermined':sum(r['decision']=='undetermined' for r in group),'reference_undetermined':sum(r['reference_decision']=='undetermined' for r in group),'incompatible':sum(not r['compatible'] for r in group)};need(s['counts'][f'{degree}/{n}/{d}']==counts,'summary mismatch')
 need(s['reference_undetermined']==sum(r['reference_decision']=='undetermined' for r in rows) and s['incompatible']==sum(not r['compatible'] for r in rows),'totals');print(json.dumps({'ok':True,'source_certificates':len(parent),'records':len(rows),'reference_undetermined':s['reference_undetermined'],'incompatible':s['incompatible'],'scope':'stored arithmetic and coverage only, not independent proof'}));return 0
if __name__=='__main__':
 try:sys.exit(main())
 except (OSError,ValueError,KeyError,TypeError,ZeroDivisionError) as e:print(json.dumps({'ok':False,'error':str(e)}));sys.exit(1)
