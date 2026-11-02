"""Independent readback of matched-grid method comparison records."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import argparse,hashlib,json,sys
MESH=(16,32,64);METHODS=('global','spatial','signed_quadratic','signed_cubic');TOL=(F(1,1000),F(1,100000))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(x,msg):
 if not x:raise ValueError(msg)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--training',required=True,type=Path);ap.add_argument('--certificates',required=True,type=Path);ap.add_argument('--directory',required=True,type=Path);a=ap.parse_args();fixture=a.training/'frozen_networks.json';h=sha(fixture);ts=json.loads((a.training/'summary.json').read_text());cs=json.loads((a.certificates/'summary.json').read_text());s=json.loads((a.directory/'summary.json').read_text());need(h==ts['fixture_sha256']==cs['fixture_sha256']==s['fixture_sha256'],'fixture hash');need(sha(a.certificates/'records.json')==cs['records_sha256']==s['input_records_sha256'] and sha(a.certificates/'references.json')==cs['references_sha256']==s['references_sha256'],'prior hash');need(sha(a.directory/'records.json')==s['records_sha256'] and sha(Path(__file__).with_name('151_new_pinn_matched_baselines_20261102.py'))==s['script_sha256'],'result/script hash');need(s['meshes']==list(MESH) and s['methods']==list(METHODS) and s['tolerances']==[str(t) for t in TOL],'config')
 nets=json.loads(fixture.read_text())['networks'];bylabel={n['label']:n for n in nets};need(len(bylabel)==len(nets)==12,'networks');refrows=json.loads((a.certificates/'references.json').read_text());refs={r['network']:r for r in refrows};need(len(refs)==len(refrows)==12 and set(refs)==set(bylabel),'references');prior={(r['network'],r['degree'],r['intervals']):r for r in json.loads((a.certificates/'records.json').read_text())};rows=json.loads((a.directory/'records.json').read_text());keys=[(r['network'],r['intervals'],r['method']) for r in rows];expected=set(product(bylabel,MESH,METHODS));need(len(rows)==len(set(keys))==len(expected)==s['records']==s['expected']==144 and set(keys)==expected,'manifest')
 for r in rows:
  net=bylabel[r['network']];q=refs[r['network']];lo,hi,rl,ru=map(F,(r['lower'],r['upper'],r['reference_lower'],r['reference_upper']));need((r['source'],r['depth'],r['seed'])==(net['source_name'],net['depth'],net['seed']),'network metadata');need(0<=lo<=hi and 0<=rl<=ru and rl==F(q['lower']) and ru==F(q['upper']) and r['reference_converged']==q['converged'] is True,'bounds/reference');need(r['reference_overlap']==(lo<=ru and rl<=hi) and r['reference_contained']==(lo<=rl and ru<=hi),'reference flags');need(r['accepted']=={str(t):hi<=t for t in TOL},'acceptance');m=r['method'];n=r['intervals']
  if m=='signed_cubic':need(r['quadratic_cache_reused'] is None and r['residual_calls']==3*n+1 and r['derivative_calls']==n,'cubic counts')
  else:need(r['quadratic_cache_reused'] is True and r['residual_calls']==2*n+1 and r['derivative_calls']==n,'quadratic counts')
  if m.startswith('signed_'):
   old=prior[(r['network'],2 if m=='signed_quadratic' else 3,n)];need(r['matches_prior']==(r['lower']==old['lower'] and r['upper']==old['upper']) and r['matches_prior'],'prior equality')
  else:need(lo==0 and r['matches_prior'] is None,'baseline lower')
 for n in MESH:
  for m in METHODS:
   g=[r for r in rows if r['intervals']==n and r['method']==m];summary={'accepted':{str(t):sum(r['accepted'][str(t)] for r in g) for t in TOL},'reference_contained':sum(r['reference_contained'] for r in g)};need(s['counts'][f'{n}/{m}']==summary,'counts')
 need(s['reference_overlap']==sum(r['reference_overlap'] for r in rows)==144 and s['reference_contained']==sum(r['reference_contained'] for r in rows) and s['shared_cache_passed']==sum(r['quadratic_cache_reused'] is True for r in rows)==108 and s['prior_matched']==sum(r['matches_prior'] is True for r in rows)==72 and not s['failures'],'summary totals');print(json.dumps({'ok':True,'records':len(rows),'references_contained':s['reference_contained'],'shared_cache_checks':s['shared_cache_passed'],'matching_prior_certificates':s['prior_matched'],'scope':'readback arithmetic and provenance, not proof of implementation'}));return 0
if __name__=='__main__':
 try:sys.exit(main())
 except (OSError,ValueError,KeyError,TypeError,ZeroDivisionError) as e:print(json.dumps({'ok':False,'error':str(e)}));sys.exit(1)
