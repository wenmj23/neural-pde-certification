"""Check hashes, exhaustive requests, budget stops, and paired medians."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from statistics import median
import argparse,hashlib,json,math,sys
MESH=(16,32,64,128,256)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(v,msg):
 if not v:raise ValueError(msg)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--training',required=True,type=Path);ap.add_argument('--certificates',required=True,type=Path);ap.add_argument('--timing',required=True,type=Path);a=ap.parse_args();t,c,p=a.training,a.certificates,a.timing;s=json.loads((p/'summary.json').read_text());ts=json.loads((t/'summary.json').read_text());cs=json.loads((c/'summary.json').read_text());fixture=t/'frozen_networks.json';need(sha(fixture)==ts['fixture_sha256']==cs['fixture_sha256']==s['fixture_sha256'],'fixture');need(sha(c/'records.json')==cs['records_sha256']==s['certificates_sha256'],'cert hash');need(sha(p/'records.json')==s['records_sha256'] and sha(p/'pairs.json')==s['pairs_sha256'],'record hash');need(sha(Path(__file__).with_name('161_sine_time_to_goal_20261102.py'))==s['script_sha256'],'script hash');nets=json.loads(fixture.read_text())['networks'];labels={n['label'] for n in nets};nmap={n['label']:n for n in nets};need(len(labels)==len(nets)==12,'labels');prior={(r['network'],r['degree'],r['intervals']):r for r in json.loads((c/'records.json').read_text())};rows=json.loads((p/'records.json').read_text());keys=[(r['network'],r['degree'],r['repeat']) for r in rows];expected=set(product(labels,(2,3),range(3)));need(len(rows)==len(set(keys))==len(expected)==s['records']==s['expected']==72 and set(keys)==expected,'records manifest');need({r['sequence'] for r in rows}==set(range(72)),'sequence')
 for r in rows:
  trace=r['trace'];need(r['frequency']==nmap[r['network']]['frequency'] and bool(trace) and len(trace)<=len(MESH) and [x['intervals'] for x in trace]==list(MESH[:len(trace)]),'mesh sequence');need(r['final_intervals']==trace[-1]['intervals'] and r['status']==('goal_met' if trace[-1]['goal'] else 'budget_exhausted') and all(not x['goal'] for x in trace[:-1]),'status');need(r['status']=='goal_met' or len(trace)==len(MESH),'budget truncation');need(type(r['wall_ns']) is int and r['wall_ns']>0 and type(r['cpu_ns']) is int and r['cpu_ns']>0,'timing')
  for x in trace:
   q=prior[(r['network'],r['degree'],x['intervals'])];lo,hi=F(x['lower']),F(x['upper']);goal=(lo==hi==0) or (lo>0 and hi<=F(51,50)*lo);need(x['goal']==goal and x['matches_prior']==(x['lower']==q['lower'] and x['upper']==q['upper']) and x['matches_prior'],'trace')
 pairs=json.loads((p/'pairs.json').read_text());need(len(pairs)==12 and {x['network'] for x in pairs}==labels,'pairs')
 for x in pairs:
  med={d:median(r['wall_ns'] for r in rows if r['network']==x['network'] and r['degree']==d) for d in (2,3)};need(x['frequency']==nmap[x['network']]['frequency'] and x['degree2_median_ns']==med[2] and x['degree3_median_ns']==med[3],'medians')
  for d in (2,3):
   subset=[r for r in rows if r['network']==x['network'] and r['degree']==d];need({r['final_intervals'] for r in subset}=={x[f'degree{d}_final_intervals']} and {r['status'] for r in subset}=={x[f'degree{d}_status']},'pair status')
  both=x['degree2_status']==x['degree3_status']=='goal_met';need((math.isclose(x['cubic_over_quadratic'],med[3]/med[2],rel_tol=1e-15) if both else x['cubic_over_quadratic'] is None),'censored ratio')
 need(s['goal_met']==sum(r['status']=='goal_met' for r in rows) and s['budget_exhausted']==sum(r['status']=='budget_exhausted' for r in rows) and s['all_trace_match']==sum(all(x['matches_prior'] for x in r['trace']) for r in rows)==72 and not s['failures'],'summary');ratios=[x['cubic_over_quadratic'] for x in pairs if x['cubic_over_quadratic'] is not None];print(json.dumps({'ok':True,'records':len(rows),'goal_met':s['goal_met'],'budget_exhausted':s['budget_exhausted'],'successful_pairs':len(ratios),'censored_pairs':len(pairs)-len(ratios),'ratio_median':median(ratios) if ratios else None,'ratio_min':min(ratios) if ratios else None,'ratio_max':max(ratios) if ratios else None,'cubic_faster':sum(x<1 for x in ratios),'scope':'readback; success-only pairs exclude budget failures and are selection-biased'}));return 0
if __name__=='__main__':
 try:sys.exit(main())
 except (OSError,ValueError,KeyError,TypeError,ZeroDivisionError) as e:print(json.dumps({'ok':False,'error':str(e)}));sys.exit(1)
