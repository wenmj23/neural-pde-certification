"""Read back the point guard output and compare to trusted input/interval run."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,sys
POWERS={'quadratic_12x2':2,'cubic_20x3':3,'quintic_42x5':5}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def key(r):return (r['fixture_sha256'],r['network'],r['source'],r['cell_index'],r['order'])
def need(b,msg):
 if not b:raise ValueError(msg)
def main():
 p=argparse.ArgumentParser();p.add_argument('--input',required=True,type=Path);p.add_argument('--points',required=True,type=Path);p.add_argument('--coarse',required=True,type=Path);p.add_argument('--fixtures',required=True,nargs='+',type=Path);a=p.parse_args()
 si=json.loads((a.input/'summary.json').read_text());sp=json.loads((a.points/'summary.json').read_text());sc=json.loads((a.coarse/'summary.json').read_text())
 for folder,s in ((a.input,si),(a.points,sp),(a.coarse,sc)):need(sha(folder/'records.json')==s['records_sha256'],'hash mismatch')
 need(sp['input_sha256']==si['records_sha256'] and sc['input_sha256']==si['records_sha256'],'input hashes');need(sha(Path(__file__).with_name('132_omission_point_guard_20261102.py'))==sp['script_sha256'],'script hash')
 labels={(sha(f),n['label']) for f in a.fixtures for n in json.loads(f.read_bytes())['networks']};expected={(h,n,s,i,j) for h,n in labels for s,pow in POWERS.items() for i in range(4) for j in (0,2,3,4) if j<=pow}
 base={key(r):r for r in json.loads((a.input/'records.json').read_text()) if r['order']<=POWERS[r['source']]};coarse={key(r):r for r in json.loads((a.coarse/'records.json').read_text()) if r['shrink']==1};rows=json.loads((a.points/'records.json').read_text());actual={key(r):r for r in rows};need(len(rows)==len(actual)==len(expected)==sp['expected']==sp['records']==864 and set(actual)==expected==set(base)==set(coarse),'manifest mismatch')
 for k,r in actual.items():
  cl,cu,bl,bu,ref=map(F,(r['correct_lower'],r['correct_upper'],r['mutated_lower'],r['mutated_upper'],r['midpoint_reference']));need(cl<=cu and bl<=bu and ref==F(base[k]['mp70']),'bounds/reference');need(r['correct_contains']==(cl<=ref<=cu) and r['correct_contains'],'correct flag');need(r['mutation_detected']==(not bl<=ref<=bu),'mutation flag')
 need(sp['correct_contains']==sum(r['correct_contains'] for r in rows) and sp['mutation_detected']==sum(r['mutation_detected'] for r in rows) and not sp['failures'],'summary');added=sum(not coarse[k]['mutation_detected'] and r['mutation_detected'] for k,r in actual.items());print(json.dumps({'ok':True,'cases':len(rows),'point_detected':sp['mutation_detected'],'coarse_detected':sum(r['mutation_detected'] for r in coarse.values()),'additional_point_detected':added,'scope':'readback consistency, not independent numerical proof'}));return 0
if __name__=='__main__':
 try:sys.exit(main())
 except (OSError,ValueError,TypeError,KeyError) as e:print(json.dumps({'ok':False,'error':str(e)}));sys.exit(1)
