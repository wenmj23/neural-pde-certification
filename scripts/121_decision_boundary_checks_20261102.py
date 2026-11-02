"""Check public classifier at exact certificate endpoints, not PDE truth."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,os,sys
from certification_pipeline import classify

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2)+'\n').encode();tmp=p.with_suffix('.tmp')
 with tmp.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 raw=a.input.read_bytes();src=json.loads(raw);rows=[]
 for index,x in enumerate(src):
  l,u=F(x['lower']),F(x['upper']);
  if not 0<=l<u:raise ValueError('expected nondegenerate bound')
  d=(u-l)/4;cases=[('at_lower',l,'undetermined'),('midpoint',(l+u)/2,'undetermined'),('at_upper',u,'within_tolerance'),('above_upper',u+d,'within_tolerance')]
  if l>0:cases.append(('below_lower',l/2,'exceeds_tolerance'))
  for kind,tau,expected in cases:
   out=classify(l,u,{'error_tolerance':str(tau),'gap_ratio':'1.02','goal':'decision'});actual=out['accuracy_decision']
   rows.append({'input_index':index,'kind':kind,'lower':str(l),'upper':str(u),'tolerance':str(tau),'expected':expected,'actual':actual,'goal_met':out['goal_met'],'passed':actual==expected and out['goal_met']==(expected!='undetermined')})
 for l,u,tau,expected in [(F(0),F(0),F(0),'within_tolerance'),(F(1,10),F(1,10),F(1,10),'within_tolerance'),(F(0),F(1,10),F(0),'undetermined')]:
  out=classify(l,u,{'error_tolerance':str(tau),'gap_ratio':'1.02','goal':'decision'});rows.append({'input_index':'synthetic','kind':'exact_zero_or_singleton','lower':str(l),'upper':str(u),'tolerance':str(tau),'expected':expected,'actual':out['accuracy_decision'],'goal_met':out['goal_met'],'passed':out['accuracy_decision']==expected and out['goal_met']==(expected!='undetermined')})
 atom(a.out/'records.json',rows);s={'input_sha256':hashlib.sha256(raw).hexdigest(),'source_certificates':len(src),'records':len(rows),'passed':sum(r['passed'] for r in rows),'records_sha256':sha(a.out/'records.json'),'script_sha256':sha(Path(__file__)),'failures':[r for r in rows if not r['passed']],'scope':'tests classifier exact comparisons using already computed bounds; not independent numerical PDE validation'};atom(a.out/'summary.json',s);print(json.dumps(s));return int(bool(s['failures']))
if __name__=='__main__':sys.exit(main())
