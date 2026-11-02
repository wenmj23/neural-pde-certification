"""Fixed-grid certificates for newly trained sine PINNs and later reference.

Uses mathematical pi in the equation; the training configuration records
that torch used its binary64 pi approximation. Reference does not drive mesh.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import argparse,hashlib,importlib.util,json,os,sys
from arb_certification import exact_fraction
from arb_cubic_reconstruction import SixthOrderDeep,CubicVerifier
from arb_shaped_remainder import ShapedVerifier
ROOT=Path(__file__).resolve().parents[1];FREQ=(1,3,5);DEPTHS=(1,2);SEEDS=(501,502);DEGREES=(2,3);MESH=(16,32,64,128,256);TOL=(F(1,1000),F(1,100),F(1,10))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2,allow_nan=False)+'\n').encode();tmp=p.with_suffix('.tmp')
 with tmp.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def decision(l,u,t):return 'within' if u<=t else ('exceeds' if l>t else 'undetermined')
def fr(a):return F(int(a['numerator']),int(a['denominator']))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--training',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False);fixture=a.training/'frozen_networks.json';h=sha(fixture);training=json.loads((a.training/'summary.json').read_text())
 if h!=training['fixture_sha256'] or training['trained']!=training['expected'] or training['failed']:raise ValueError('training incomplete or fixture changed')
 nets=json.loads(fixture.read_bytes())['networks'];expected=set(product(FREQ,DEPTHS,SEEDS));actual={(n['frequency'],n['depth'],n['seed']) for n in nets}
 if len(nets)!=len(expected) or len(actual)!=len(expected) or actual!=expected:raise ValueError('network manifest')
 path=ROOT/'scripts/24_depth_continuous_reference.py';spec=importlib.util.spec_from_file_location('reference24',path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
 cfg={'fixture_sha256':h,'network_count':len(nets),'frequencies':FREQ,'meshes':MESH,'degrees':DEGREES,'tolerances':[str(t) for t in TOL],'gap_target':'51/50','reference':'after both degrees and all meshes; mpmath interval relative 1e-7, absolute 1e-12, cap 10000','script_sha256':sha(Path(__file__)),'reference_script_sha256':sha(path),'scope':'finite trained sine PINNs; conditional on smooth model and interval backends'};atom(a.out/'locked_config.json',cfg)
 rows=[];refs=[]
 for net in nets:
  pair=[]
  for degree in DEGREES:
   p=SixthOrderDeep(fixture,net['label'],h);v=ShapedVerifier(p,'moment') if degree==2 else CubicVerifier(p,'moment')
   for n in MESH:
    v.uniform(n);r=v.certify();lo,hi=exact_fraction(r['lower']),exact_fraction(r['upper']);zl,zu,du=(exact_fraction(r[x]) for x in ('z_lower','z_upper','remainder_upper'))
    if not 0<=lo<=hi:raise ValueError('invalid interval')
    pair.append({'network':net['label'],'frequency':net['frequency'],'depth':net['depth'],'seed':net['seed'],'degree':degree,'intervals':n,'lower':str(lo),'upper':str(hi),'z_lower':str(zl),'z_upper':str(zu),'remainder_upper':str(du),'gap_1_02':bool((lo==hi==0) or (lo>0 and hi<=F(51,50)*lo)),'decisions':{str(t):decision(lo,hi,t) for t in TOL}})
  ref=mod.MPDeepReference(net).maximum(relative=F(1,10**7),absolute=F(1,10**12),cap=10000);rl,ru=fr(ref['lower']),fr(ref['upper']);q={'network':net['label'],'lower':str(rl),'upper':str(ru),'converged':bool(ref['converged']),'splits':ref['splits'],'decisions':{str(t):decision(rl,ru,t) for t in TOL}};refs.append(q)
  for r in pair:
   lo,hi=F(r['lower']),F(r['upper']);r.update(reference_converged=q['converged'],reference_overlap=lo<=ru and rl<=hi,reference_contained=lo<=rl and ru<=hi,decision_compatible=all(r['decisions'][str(t)]=='undetermined' or q['decisions'][str(t)]=='undetermined' or r['decisions'][str(t)]==q['decisions'][str(t)] for t in TOL));rows.append(r)
  print(net['label'],'reference converged',q['converged'],flush=True)
 atom(a.out/'records.json',rows);atom(a.out/'references.json',refs);counts={}
 for k in FREQ:
  for d in DEGREES:
   for n in MESH:
    g=[r for r in rows if r['frequency']==k and r['degree']==d and r['intervals']==n];counts[f'{k}/{d}/{n}']={'gap_met':sum(r['gap_1_02'] for r in g),'overlap':sum(r['reference_overlap'] for r in g),'contained':sum(r['reference_contained'] for r in g),'decisions':{str(t):{status:sum(r['decisions'][str(t)]==status for r in g) for status in ('within','exceeds','undetermined')} for t in TOL}}
 s={'fixture_sha256':h,'expected':len(nets)*len(DEGREES)*len(MESH),'records':len(rows),'references':len(refs),'reference_converged':sum(r['converged'] for r in refs),'overlap':sum(r['reference_overlap'] for r in rows),'contained':sum(r['reference_contained'] for r in rows),'compatible':sum(r['decision_compatible'] for r in rows),'counts':counts,'records_sha256':sha(a.out/'records.json'),'references_sha256':sha(a.out/'references.json'),'config_sha256':sha(a.out/'locked_config.json'),'failures':[{'network':r['network'],'degree':r['degree'],'mesh':r['intervals'],'reference_converged':r['reference_converged'],'overlap':r['reference_overlap']} for r in rows if not r['reference_converged'] or not r['reference_overlap'] or not r['decision_compatible']]}
 if len(rows)!=s['expected']:raise ValueError('record count')
 atom(a.out/'summary.json',s);print(json.dumps(s),flush=True);return int(bool(s['failures']))
if __name__=='__main__':sys.exit(main())
