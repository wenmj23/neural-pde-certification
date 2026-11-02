"""Read back pre-specified full-domain and sine stress result records."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,sys

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(v,msg):
 if not v:raise ValueError(msg)
def verify(p,kind,fixtures):
 s=json.loads((p/'summary.json').read_text());rows=json.loads((p/'records.json').read_text());need(isinstance(rows,list),'records not list')
 need(digest(p/'records.json')==s['records_sha256' if kind=='full' else 'sha256'],'hash mismatch')
 need(len(rows)==s['records']==s['expected'],'count mismatch')
 if kind=='full':
  labels={(digest(f),n['label']) for f in fixtures for n in json.loads(f.read_bytes())['networks']};keys=[(r['fixture_sha256'],r['network'],r['cell_index'],r['order']) for r in rows]
  expected={(sha,label,i,o) for sha,label in labels for i in range(32) for o in (0,2,3,4)}
  need(len(keys)==len(expected) and set(keys)==expected,'fixture/coverage/duplicate mismatch')
  for r in rows:
   lo,hi,il,iu=map(F,(r['arb_lower'],r['arb_upper'],r['mpiv_lower'],r['mpiv_upper']))
   need(lo<=hi and il<=iu,'reversed interval');need(r['interval_overlap']==(lo<=iu and il<=hi),'incorrect overlap flag')
  for key in ('interval_overlap','midpoint_stable','midpoint_in_arb','midpoint_in_mpiv'):
   need(all(type(r[key]) is bool for r in rows),'flag not boolean');count=sum(r[key] for r in rows);need(s[key]==count==len(rows),'failed or incorrect '+key)
 else:
  frequencies,meshes=([1,4,16,32],[16,32,64]) if kind=='sine' else ([16,32],[128,256,512])
  keys=[(r['frequency'],r['intervals'],r['degree']) for r in rows];expected={(f,n,d) for f in frequencies for n in meshes for d in (2,3)}
  need(len(keys)==len(expected) and set(keys)==expected,'sine key coverage mismatch')
  for r in rows:
   lo,hi,truth=map(F,(r['lower'],r['upper'],r['known_exact_error']));need(truth==F(1,10000) and 0<=lo<=hi,'invalid sine bounds')
   need(r['contains_truth']==(lo<=truth<=hi) and r['contains_truth'],'truth coverage error')
   need(r['gap_goal_1_02']==(lo>0 and hi<=F(51,50)*lo),'incorrect ratio flag')
  need(s['contains_truth']==len(rows) and s['goal_met']==sum(r['gap_goal_1_02'] for r in rows),'summary mismatch')
 need(not s['failures'],'reported failures')
 return {'ok':True,'kind':kind,'records':len(rows),'sha256':digest(p/'records.json'),'scope':'record integrity and stored consistency, not correctness of numerical calculation'}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('kind',choices=['full','sine','refinement']);ap.add_argument('directory',type=Path);ap.add_argument('--fixtures',nargs='*',type=Path,default=[]);a=ap.parse_args()
 try:print(json.dumps(verify(a.directory,a.kind,a.fixtures)));return 0
 except (OSError,ValueError,KeyError,TypeError,ZeroDivisionError) as e:print(json.dumps({'ok':False,'error':str(e)}));return 1
if __name__=='__main__':sys.exit(main())
