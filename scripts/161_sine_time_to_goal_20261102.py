"""Exploratory sine PINN time to 1.02 enclosure goal, preserving budget failures.

Same network/degree schedule, 1 warmup + 3 fresh instances per method; timings
include prior attempted meshes and object/fixture setup, not Python startup.
"""
from pathlib import Path
from fractions import Fraction as F
from statistics import median
import argparse,hashlib,json,os,random,sys,time
from arb_certification import exact_fraction
from arb_cubic_reconstruction import SixthOrderDeep,CubicVerifier
from arb_shaped_remainder import ShapedVerifier
MESH=(16,32,64,128,256);DEGREES=(2,3);REPEATS=3;SEED=20261105
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2,allow_nan=False)+'\n').encode();tmp=p.with_suffix('.tmp')
 with tmp.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--training',required=True,type=Path);ap.add_argument('--certificates',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False);fixture=a.training/'frozen_networks.json';h=sha(fixture);ts=json.loads((a.training/'summary.json').read_text());cs=json.loads((a.certificates/'summary.json').read_text());cp=a.certificates/'records.json'
 if h!=ts['fixture_sha256'] or h!=cs['fixture_sha256'] or sha(cp)!=cs['records_sha256']:raise ValueError('input hash')
 nets=json.loads(fixture.read_bytes())['networks'];prior={(r['network'],r['degree'],r['intervals']):r for r in json.loads(cp.read_text())}
 if len(prior)!=len(nets)*len(DEGREES)*len(MESH):raise ValueError('prior manifest')
 def run(net,degree):
  start=time.perf_counter_ns();cpu=time.process_time_ns();p=SixthOrderDeep(fixture,net['label'],h);v=ShapedVerifier(p,'moment') if degree==2 else CubicVerifier(p,'moment');trace=[];goal=False
  for n in MESH:
   v.uniform(n);r=v.certify();lo,hi=exact_fraction(r['lower']),exact_fraction(r['upper']);goal=(lo==hi==0) or (lo>0 and hi<=F(51,50)*lo);q=prior[(net['label'],degree,n)];trace.append({'intervals':n,'lower':str(lo),'upper':str(hi),'goal':goal,'matches_prior':str(lo)==q['lower'] and str(hi)==q['upper']})
   if goal:break
  cpu=time.process_time_ns()-cpu;wall=time.perf_counter_ns()-start;return wall,cpu,goal,trace
 for net in nets:
  for d in DEGREES:run(net,d)
 tasks=[(net,d,rep) for net in nets for d in DEGREES for rep in range(REPEATS)];random.Random(SEED).shuffle(tasks);rows=[]
 for seq,(net,d,rep) in enumerate(tasks):
  wall,cpu,goal,trace=run(net,d);rows.append({'sequence':seq,'network':net['label'],'frequency':net['frequency'],'degree':d,'repeat':rep,'wall_ns':wall,'cpu_ns':cpu,'status':'goal_met' if goal else 'budget_exhausted','final_intervals':trace[-1]['intervals'],'trace':trace})
 atom(a.out/'records.json',rows);pairs=[]
 for n in nets:
  groups={d:[r for r in rows if r['network']==n['label'] and r['degree']==d] for d in DEGREES};met={d:all(r['status']=='goal_met' for r in groups[d]) for d in DEGREES};med={d:median(r['wall_ns'] for r in groups[d]) for d in DEGREES};mesh={d:{r['final_intervals'] for r in groups[d]} for d in DEGREES}
  if any(len(mesh[d])!=1 for d in DEGREES):raise ValueError('final mesh varies across repeats')
  pairs.append({'network':n['label'],'frequency':n['frequency'],'degree2_median_ns':med[2],'degree3_median_ns':med[3],'degree2_status':'goal_met' if met[2] else 'budget_exhausted','degree3_status':'goal_met' if met[3] else 'budget_exhausted','degree2_final_intervals':next(iter(mesh[2])),'degree3_final_intervals':next(iter(mesh[3])),'cubic_over_quadratic':med[3]/med[2] if all(met.values()) else None})
 atom(a.out/'pairs.json',pairs);summary={'fixture_sha256':h,'certificates_sha256':cs['records_sha256'],'seed':SEED,'meshes':MESH,'degrees':DEGREES,'repeats':REPEATS,'expected':len(nets)*len(DEGREES)*REPEATS,'records':len(rows),'goal_met':sum(r['status']=='goal_met' for r in rows),'budget_exhausted':sum(r['status']=='budget_exhausted' for r in rows),'all_trace_match':sum(all(x['matches_prior'] for x in r['trace']) for r in rows),'records_sha256':sha(a.out/'records.json'),'pairs_sha256':sha(a.out/'pairs.json'),'script_sha256':sha(Path(__file__)),'failures':[{'network':r['network'],'degree':r['degree']} for r in rows if not all(x['matches_prior'] for x in r['trace']) or r['wall_ns']<=0 or r['cpu_ns']<=0],'scope':'post-hoc one-host cumulative goal timing; budget failures retained and excluded from success-only time ratios'}
 if len(rows)!=summary['expected']:raise ValueError('record count')
 atom(a.out/'summary.json',summary);print(json.dumps(summary),flush=True);return int(bool(summary['failures']))
if __name__=='__main__':sys.exit(main())
