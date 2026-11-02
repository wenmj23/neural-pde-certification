"""Separate-process integrity and consistency check for scripts 143 and 144."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import argparse,hashlib,json,math,sys
SOURCES={'quadratic_12x2':[0,0,12],'cubic_20x3':[0,0,0,20],'quintic_42x5':[0,0,0,0,0,42]};MESH=(16,32,64,128,256);DEGREES=(2,3);TOL=(F(1,1000),F(1,100000))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(v,msg):
 if not v:raise ValueError(msg)
def decision(l,u,t):return 'within' if u<=t else ('exceeds' if l>t else 'undetermined')
def main():
 p=argparse.ArgumentParser();p.add_argument('--training',type=Path,required=True);p.add_argument('--certificates',type=Path,required=True);a=p.parse_args();t=a.training;c=a.certificates
 ts=json.loads((t/'summary.json').read_text());tc=json.loads((t/'locked_config.json').read_text());logs=json.loads((t/'training_log.json').read_text());fixture=json.loads((t/'frozen_networks.json').read_text());h=sha(t/'frozen_networks.json');need(h==ts['fixture_sha256'] and sha(t/'training_log.json')==ts['training_log_sha256'] and sha(t/'locked_config.json')==ts['config_sha256'],'training hash');need(sha(Path(__file__).with_name('143_train_new_forcing_20261102.py'))==tc['script_sha256'],'train script');nets=fixture['networks'];keys={(n['source_name'],n['depth'],n['seed']) for n in nets};expected=set(product(SOURCES,(1,2),(401,402)));need(len(nets)==len(keys)==ts['trained']==ts['expected']==12 and keys==expected and not ts['failed'],'training manifest');need(tc['sources']==SOURCES and fixture['config']==tc,'training config')
 for n in nets:need(n['source_coefficients']==SOURCES[n['source_name']] and len(n['layers'])==n['depth']+1,'training source/depth')
 labels={n['label'] for n in nets};need(len(labels)==12 and len(logs)==ts['log_records']==12 and {x['label'] for x in logs}==labels and all(x['status']=='trained' and math.isfinite(x['final_collocation_loss']) for x in logs),'training logs')
 s=json.loads((c/'summary.json').read_text());cfg=json.loads((c/'locked_config.json').read_text());need(h==s['fixture_sha256']==cfg['fixture_sha256'],'fixture link');need(sha(c/'locked_config.json')==s['config_sha256'] and sha(c/'records.json')==s['records_sha256'] and sha(c/'references.json')==s['references_sha256'],'cert hash');need(sha(Path(__file__).with_name('144_new_forcing_certify_20261102.py'))==cfg['script_sha256'],'cert script');need(cfg['meshes']==list(MESH) and cfg['degrees']==list(DEGREES) and cfg['tolerances']==[str(t) for t in TOL] and cfg['sources']==SOURCES,'cert config')
 refs=json.loads((c/'references.json').read_text());ref={r['network']:r for r in refs};need(len(refs)==len(ref)==s['references']==12 and set(ref)==labels,'reference manifest');rows=json.loads((c/'records.json').read_text());actual=[(r['network'],r['degree'],r['intervals']) for r in rows];expected_rows=set(product(labels,DEGREES,MESH));need(len(rows)==len(set(actual))==s['records']==s['expected']==120 and set(actual)==expected_rows,'certificate manifest');nmap={n['label']:n for n in nets}
 for r in rows:
  n=nmap[r['network']];q=ref[r['network']];need((r['source'],r['depth'],r['seed'])==(n['source_name'],n['depth'],n['seed']),'row network metadata');l,u,zl,zu,du,rl,ru=map(F,(r['lower'],r['upper'],r['z_lower'],r['z_upper'],r['remainder_upper'],q['lower'],q['upper']));need(0<=l<=u and 0<=zl<=zu and du>=0 and l<=max(F(0),zl-du) and u>=zu+du and 0<=rl<=ru,'invalid arithmetic');need(r['gap_1_02']==((l==u==0) or (l>0 and u<=F(51,50)*l)),'gap');need(r['reference_converged']==q['converged'] and r['reference_overlap']==(l<=ru and rl<=u) and r['reference_contained']==(l<=rl and ru<=u),'reference flags');need(set(r['decisions'])=={str(t) for t in TOL} and set(q['decisions'])=={str(t) for t in TOL},'decision keys')
  for tol in TOL:need(r['decisions'][str(tol)]==decision(l,u,tol) and q['decisions'][str(tol)]==decision(rl,ru,tol),'decision')
  compat=all(r['decisions'][str(tol)]=='undetermined' or q['decisions'][str(tol)]=='undetermined' or r['decisions'][str(tol)]==q['decisions'][str(tol)] for tol in TOL);need(r['decision_compatible']==compat,'compat flag')
 for d in DEGREES:
  for mesh in MESH:
   g=[r for r in rows if r['degree']==d and r['intervals']==mesh];counts={'gap_met':sum(r['gap_1_02'] for r in g),'overlap':sum(r['reference_overlap'] for r in g),'contained':sum(r['reference_contained'] for r in g),'decisions':{str(t):{k:sum(r['decisions'][str(t)]==k for r in g) for k in ('within','exceeds','undetermined')} for t in TOL}};need(s['counts'][f'{d}/{mesh}']==counts,'summary group')
 need(s['reference_converged']==sum(r['converged'] for r in refs) and s['overlap']==sum(r['reference_overlap'] for r in rows) and s['contained']==sum(r['reference_contained'] for r in rows) and s['compatible']==sum(r['decision_compatible'] for r in rows),'summary totals');need(not s['failures'],'reported failures');print(json.dumps({'ok':True,'networks':len(nets),'certificates':len(rows),'references':len(refs),'reference_converged':s['reference_converged'],'overlap':s['overlap'],'contained':s['contained'],'compatible':s['compatible'],'scope':'storage consistency and exact saved arithmetic; not a proof of numerical computation'}));return 0
if __name__=='__main__':
 try:sys.exit(main())
 except (OSError,KeyError,TypeError,ValueError,ZeroDivisionError) as e:print(json.dumps({'ok':False,'error':str(e)}));sys.exit(1)
