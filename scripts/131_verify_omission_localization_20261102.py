"""Read-back validation for fault-injection experiments 127-130."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,sys

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(b,msg):
 if not b:raise ValueError(msg)
def key(r):return (r['fixture_sha256'],r['network'],r['source'],r['cell_index'],r['order'])
def load(p,script):
 s=json.loads((p/'summary.json').read_text());need(sha(p/'records.json')==s['records_sha256'],'record hash '+str(p));need(sha(Path(__file__).parent/script)==s['script_sha256'],'script hash '+str(p));rows=json.loads((p/'records.json').read_text());need(len(rows)==s['records']==s['expected'],'count '+str(p));need(not s['failures'],'failures '+str(p));return s,rows
def validate(rows,summary,keys,refmap,parameter,values):
 seen=set()
 for r in rows:
  k=key(r);v=r[parameter];need(k in keys and v in values and (k,v) not in seen,'unexpected or duplicate key');seen.add((k,v));cl,cu,bl,bu,ref=map(F,(r['correct_lower'],r['correct_upper'],r['mutated_lower'],r['mutated_upper'],r['midpoint_reference']));need(cl<=cu and bl<=bu,'reversed bounds');need(ref==refmap[k],'reference mismatch');need(r['correct_contains']==(cl<=ref<=cu) and r['correct_contains'],'correct flag');need(r['mutation_detected']==(not bl<=ref<=bu),'mutation flag')
 need(seen=={(k,v) for k in keys for v in values},'coverage mismatch')
 for v in values:
  subset=[r for r in rows if r[parameter]==v];got=summary['by_'+('bits' if parameter=='bits' else 'shrink')][str(v)];need(got['total']==len(subset) and got['correct_contains']==sum(r['correct_contains'] for r in subset) and got['mutation_detected']==sum(r['mutation_detected'] for r in subset),'summary mismatch')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--fixtures',type=Path,nargs='+',required=True);a=ap.parse_args();p=a.input
 base=json.loads((p/'records.json').read_text());bs=json.loads((p/'summary.json').read_text());need(sha(p/'records.json')==bs['records_sha256'],'base hash');labels={(sha(f),n['label']) for f in a.fixtures for n in json.loads(f.read_bytes())['networks']};powers={'quadratic_12x2':2,'cubic_20x3':3,'quintic_42x5':5};base={key(r):r for r in base if r['order']<=powers[r['source']]};expected={(h,n,src,i,j) for h,n in labels for src,pow in powers.items() for i in range(4) for j in (0,2,3,4) if j<=pow};need(len(base)==len(expected)==864 and set(base)==expected,'base keys');ref={k:F(r['mp70']) for k,r in base.items()}
 s1,r1=load(a.root/'coarse','127_mutated_forcing_localization_20261102.py');need(s1['input_sha256']==bs['records_sha256'],'coarse input');validate(r1,s1,expected,ref,'shrink',[1,2,4,8,16]);keys16={key(r) for r in r1 if r['shrink']==16 and not r['mutation_detected']}
 s2,r2=load(a.root/'finer','128_mutation_finer_followup_20261102.py');need(s2['input_sha256']==s1['records_sha256'] and s2['selected_count']==len(keys16),'finer input');validate(r2,s2,keys16,ref,'shrink',[32,64,128,256,1024,4096]);keys4096={key(r) for r in r2 if r['shrink']==4096 and not r['mutation_detected']}
 s3,r3=load(a.root/'edge','129_mutation_edge_followup_20261102.py');need(s3['input_sha256']==s2['records_sha256'] and s3['selected_count']==len(keys4096),'edge input');validate(r3,s3,keys4096,ref,'shrink',[16384,65536,262144,1048576,4194304]);s4,r4=load(a.root/'point','130_point_precision_omission_20261102.py');need(s4['input_sha256']==s2['records_sha256'] and s4['selected_count']==len(keys4096),'point input');validate(r4,s4,keys4096,ref,'bits',[103,160]);print(json.dumps({'ok':True,'base_cases':len(expected),'coarse':len(r1),'finer':len(r2),'edge':len(r3),'point':len(r4),'scope':'readback consistency; numerical calculations are not rerun'}));return 0
if __name__=='__main__':
 try:sys.exit(main())
 except (ValueError,KeyError,TypeError,OSError) as e:print(json.dumps({'ok':False,'error':str(e)}));sys.exit(1)
