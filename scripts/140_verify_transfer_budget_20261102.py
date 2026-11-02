"""Independently validate stored keys, exact decisions and input provenance."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,sys
SOURCES={'quadratic_12x2':[0,0,12],'cubic_20x3':[0,0,0,20],'quintic_42x5':[0,0,0,0,0,42]};MESH=(8,16,32,64);DEGREES=(2,3);TOL=(F(1,10),F(1,5),F(3,10))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(b,msg):
 if not b:raise ValueError(msg)
def key(r):return (r['fixture_sha256'],r['network'],r['source'])
def classify(l,u,t):return 'within' if u<=t else ('exceeds' if l>t else 'undetermined')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--fixtures',nargs='+',required=True,type=Path);ap.add_argument('--references',required=True,type=Path);ap.add_argument('--directory',required=True,type=Path);a=ap.parse_args();s=json.loads((a.directory/'summary.json').read_text());rpath=a.directory/'records.json';need(sha(rpath)==s['records_sha256'],'record hash');need(sha(Path(__file__).with_name('139_transfer_decision_budget_20261102.py'))==s['script_sha256'],'script hash')
 old=json.loads((a.references/'summary.json').read_text());refpath=a.references/'references.json';need(sha(refpath)==old['references_sha256']==s['reference_sha256'],'reference hash');refs=json.loads(refpath.read_text());refmap={key(r):r for r in refs};labels={(sha(f),n['label']) for f in a.fixtures for n in json.loads(f.read_bytes())['networks']};basekeys={(h,n,src) for h,n in labels for src in SOURCES};need(len(refs)==len(refmap)==len(basekeys) and set(refmap)==basekeys,'reference manifest');need(s['sources']==SOURCES and s['degrees']==list(DEGREES) and s['meshes']==list(MESH) and s['tolerances']==[str(t) for t in TOL],'config mismatch')
 records=json.loads(rpath.read_text());keys=[(*key(r),r['degree'],r['intervals']) for r in records];expected={(*k,d,n) for k in basekeys for d in DEGREES for n in MESH};need(len(records)==len(set(keys))==len(expected)==s['records']==s['expected'] and set(keys)==expected,'record manifest')
 refrows=json.loads((a.references/'records.json').read_text());need(sha(a.references/'records.json')==old['records_sha256'],'prior certificates hash');previous={(*key(r),r['degree']):r for r in refrows};matches=0
 for r in records:
  ref=refmap[key(r)];l,u,rl,ru=map(F,(r['lower'],r['upper'],r['reference_lower'],r['reference_upper']));need(0<=l<=u and ref['converged'] is True and rl==F(ref['lower']) and ru==F(ref['upper']),'reference mismatch');need(r['reference_overlap']==(l<=ru and rl<=u) and r['reference_overlap'],'reference overlap');need(r['gap_1_02']==(l>0 and u<=F(51,50)*l),'gap mismatch');need(set(r['decisions'])=={str(t) for t in TOL},'tolerance keys')
  for t in TOL:
   v=r['decisions'][str(t)];d,refd=classify(l,u,t),classify(rl,ru,t);compatible=d=='undetermined' or refd=='undetermined' or d==refd;need(v['decision']==d and v['reference_decision']==refd and v['compatible']==compatible and compatible,'decision mismatch')
  if r['intervals']==64:
   prev=previous[(*key(r),r['degree'])];matches+=(r['lower']==prev['lower'] and r['upper']==prev['upper'])
 for d in DEGREES:
  for n in MESH:
   group=[r for r in records if r['degree']==d and r['intervals']==n]
   for t in TOL:
    vals=[r['decisions'][str(t)] for r in group];counts={'within':sum(x['decision']=='within' for x in vals),'exceeds':sum(x['decision']=='exceeds' for x in vals),'undetermined':sum(x['decision']=='undetermined' for x in vals),'reference_undetermined':sum(x['reference_decision']=='undetermined' for x in vals),'incompatible':sum(not x['compatible'] for x in vals)};need(s['counts'][f'{d}/{n}/{t}']==counts,'summary mismatch')
 need(not s['failures'],'reported failures');print(json.dumps({'ok':True,'records':len(records),'decisions':len(records)*len(TOL),'overlap':sum(r['reference_overlap'] for r in records),'matching_prior_64_intervals':matches,'expected_prior_comparisons':len(basekeys)*len(DEGREES),'scope':'readback consistency and reference compatibility, not mathematical audit'}));return 0
if __name__=='__main__':
 try:sys.exit(main())
 except (OSError,KeyError,TypeError,ValueError,ZeroDivisionError) as e:print(json.dumps({'ok':False,'error':str(e)}));sys.exit(1)
