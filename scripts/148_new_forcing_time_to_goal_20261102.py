"""Matched enclosure-goal timing, including every attempted grid.

Single process/host; one warmup and three randomized fresh verifier runs per
network/degree. Only fixed mesh schedule, no time-to-decision conclusion.
"""
from pathlib import Path
from fractions import Fraction as F
from statistics import median
import argparse,hashlib,json,os,random,sys,time
from arb_certification import exact_fraction
from arb_cubic_reconstruction import CubicVerifier
from arb_shaped_remainder import ShapedVerifier
ROOT=Path(__file__).resolve().parents[1];MESH=(16,32,64,128,256);DEGREES=(2,3);REPEATS=3;SEED=20261103
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2,allow_nan=False)+'\n').encode();t=p.with_suffix('.tmp')
 with t.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(t,p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--training',required=True,type=Path);ap.add_argument('--certificates',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 fixture=a.training/'frozen_networks.json';h=sha(fixture);ts=json.loads((a.training/'summary.json').read_text());cs=json.loads((a.certificates/'summary.json').read_text());cp=a.certificates/'records.json'
 if h!=ts['fixture_sha256'] or h!=cs['fixture_sha256'] or sha(cp)!=cs['records_sha256']:raise ValueError('input hash')
 nets=json.loads(fixture.read_bytes())['networks'];prior={(r['network'],r['degree'],r['intervals']):r for r in json.loads(cp.read_text())}
 if len(prior)!=len(nets)*len(DEGREES)*len(MESH):raise ValueError('prior manifest')
 src=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
 if src.count(needle)!=1:raise ValueError('prototype changed')
 ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(src.replace(needle,'OUT=None'),'defs','exec'),ns)
 def run(net,degree):
  start=time.perf_counter_ns();cpu=time.process_time_ns();p=ns['PolynomialForcing'](fixture,net['label'],h,net['source_coefficients']);v=ShapedVerifier(p,'moment') if degree==2 else CubicVerifier(p,'moment');trace=[];goal=False
  for n in MESH:
   v.uniform(n);r=v.certify();l,u=exact_fraction(r['lower']),exact_fraction(r['upper']);goal=(l==u==0) or (l>0 and u<=F(51,50)*l);q=prior[(net['label'],degree,n)];trace.append({'intervals':n,'lower':str(l),'upper':str(u),'goal':goal,'matches_prior':str(l)==q['lower'] and str(u)==q['upper']})
   if goal:break
  cpu=time.process_time_ns()-cpu;elapsed=time.perf_counter_ns()-start;return elapsed,cpu,goal,trace
 for net in nets:
  for degree in DEGREES:run(net,degree)
 tasks=[(n,d,i) for n in nets for d in DEGREES for i in range(REPEATS)];random.Random(SEED).shuffle(tasks);rows=[]
 for i,(net,d,rep) in enumerate(tasks):
  wall,cpu,goal,trace=run(net,d);rows.append({'sequence':i,'network':net['label'],'degree':d,'repeat':rep,'wall_ns':wall,'cpu_ns':cpu,'goal_met':goal,'final_intervals':trace[-1]['intervals'],'trace':trace})
 atom(a.out/'records.json',rows);pairs=[]
 for net in nets:
  med={d:median(r['wall_ns'] for r in rows if r['network']==net['label'] and r['degree']==d) for d in DEGREES};grid={d:{r['final_intervals'] for r in rows if r['network']==net['label'] and r['degree']==d} for d in DEGREES}
  if any(len(grid[d])!=1 for d in DEGREES):raise ValueError('grid inconsistency')
  pairs.append({'network':net['label'],'degree2_median_ns':med[2],'degree3_median_ns':med[3],'ratio_cubic_over_quadratic':med[3]/med[2],'degree2_final_intervals':next(iter(grid[2])),'degree3_final_intervals':next(iter(grid[3]))})
 atom(a.out/'pairs.json',pairs);summary={'fixture_sha256':h,'certificates_sha256':cs['records_sha256'],'seed':SEED,'meshes':MESH,'degrees':DEGREES,'repeats':REPEATS,'target':'U/L<=1.02, or exact zero','expected':len(nets)*len(DEGREES)*REPEATS,'records':len(rows),'goal_met':sum(r['goal_met'] for r in rows),'all_trace_match':sum(all(t['matches_prior'] for t in r['trace']) for r in rows),'records_sha256':sha(a.out/'records.json'),'pairs_sha256':sha(a.out/'pairs.json'),'script_sha256':sha(Path(__file__)),'scope':'one host, one warmup, three repeats, cumulative 16-to-target schedule, excludes Python startup/training/reference; not universal','failures':[{'network':r['network'],'degree':r['degree']} for r in rows if not r['goal_met'] or not all(t['matches_prior'] for t in r['trace']) or r['wall_ns']<=0 or r['cpu_ns']<=0]}
 if len(rows)!=summary['expected']:raise ValueError('count mismatch')
 atom(a.out/'summary.json',summary);print(json.dumps(summary),flush=True);return int(bool(summary['failures']))
if __name__=='__main__':sys.exit(main())
