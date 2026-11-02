"""Fixed mesh/degree certification of locally trained networks.

Targets, meshes and degree choices are fixed before certificate/reference runs.
Reference is calculated only after all ten certificate attempts for a network.
No timing comparisons, no automatic publication of generated weights.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,os,sys
from arb_certification import exact_fraction
from arb_cubic_reconstruction import CubicVerifier
from arb_shaped_remainder import ShapedVerifier
ROOT=Path(__file__).resolve().parents[1];SOURCES={'quadratic_12x2':[0,0,12],'cubic_20x3':[0,0,0,20],'quintic_42x5':[0,0,0,0,0,42]};DEGREES=(2,3);MESH=(16,32,64,128,256);TOL=(F(1,1000),F(1,100000))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2,allow_nan=False)+'\n').encode();t=p.with_suffix('.tmp')
 with t.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(t,p)
def decision(l,u,t):return 'within' if u<=t else ('exceeds' if l>t else 'undetermined')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--training',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 training=json.loads((a.training/'summary.json').read_text());fixture=a.training/'frozen_networks.json';h=sha(fixture)
 if h!=training['fixture_sha256'] or training['trained']!=training['expected'] or training['failed']:raise ValueError('incomplete or modified training fixture')
 nets=json.loads(fixture.read_bytes())['networks'];expected={(s,d,seed) for s in SOURCES for d in (1,2) for seed in (401,402)};actual={(n['source_name'],n['depth'],n['seed']) for n in nets}
 if len(nets)!=len(actual)==len(expected) or actual!=expected:raise ValueError('unexpected network set')
 for n in nets:
  if n['source_coefficients']!=SOURCES[n['source_name']]:raise ValueError('source mismatch')
 src=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
 if src.count(needle)!=1:raise ValueError('prototype changed')
 ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(src.replace(needle,'OUT=None'),'defs','exec'),ns)
 config={'fixture_sha256':h,'fixture_records':len(nets),'sources':SOURCES,'meshes':MESH,'degrees':DEGREES,'tolerances':[str(t) for t in TOL],'enclosure_ratio':'51/50','reference':'mpmath interval, computed after both degree grids; relative 1e-7, absolute 1e-14, cap 5000','script_sha256':sha(Path(__file__)),'scope':'finite new PINNs; conditional on interval implementation and mathematical assumptions'};atom(a.out/'locked_config.json',config)
 rows=[];refs=[]
 for n in nets:
  pair=[]
  for degree in DEGREES:
   p=ns['PolynomialForcing'](fixture,n['label'],h,n['source_coefficients']);v=ShapedVerifier(p,'moment') if degree==2 else CubicVerifier(p,'moment')
   for mesh in MESH:
    v.uniform(mesh);r=v.certify();l,u=exact_fraction(r['lower']),exact_fraction(r['upper']);zl,zu,du=(exact_fraction(r[k]) for k in ('z_lower','z_upper','remainder_upper'))
    if not 0<=l<=u:raise ValueError('invalid interval')
    row={'network':n['label'],'source':n['source_name'],'depth':n['depth'],'seed':n['seed'],'degree':degree,'intervals':mesh,'lower':str(l),'upper':str(u),'z_lower':str(zl),'z_upper':str(zu),'remainder_upper':str(du),'gap_1_02':bool((l==u==0) or (l>0 and u<=F(51,50)*l)),'decisions':{str(t):decision(l,u,t) for t in TOL}};pair.append(row)
  ref=ns['PolynomialReference'](n,n['source_coefficients']).maximum(relative=F(1,10**7),absolute=F(1,10**14),cap=5000);rl,ru=ns['fr'](ref['lower']),ns['fr'](ref['upper']);rr={'network':n['label'],'lower':str(rl),'upper':str(ru),'converged':bool(ref['converged']),'splits':ref['splits'],'decisions':{str(t):decision(rl,ru,t) for t in TOL}};refs.append(rr)
  for row in pair:
   l,u=F(row['lower']),F(row['upper']);row.update(reference_converged=rr['converged'],reference_overlap=l<=ru and rl<=u,reference_contained=l<=rl and ru<=u,decision_compatible=all(row['decisions'][str(t)]=='undetermined' or rr['decisions'][str(t)]=='undetermined' or row['decisions'][str(t)]==rr['decisions'][str(t)] for t in TOL));rows.append(row)
  print(n['label'],'completed ref converged',rr['converged'],flush=True)
 atom(a.out/'records.json',rows);atom(a.out/'references.json',refs);counts={}
 for degree in DEGREES:
  for mesh in MESH:
   g=[r for r in rows if r['degree']==degree and r['intervals']==mesh];counts[f'{degree}/{mesh}']={'gap_met':sum(r['gap_1_02'] for r in g),'overlap':sum(r['reference_overlap'] for r in g),'contained':sum(r['reference_contained'] for r in g),'decisions':{str(t):{k:sum(r['decisions'][str(t)]==k for r in g) for k in ('within','exceeds','undetermined')} for t in TOL}}
 summary={'fixture_sha256':h,'expected':len(nets)*len(DEGREES)*len(MESH),'records':len(rows),'references':len(refs),'reference_converged':sum(r['converged'] for r in refs),'overlap':sum(r['reference_overlap'] for r in rows),'contained':sum(r['reference_contained'] for r in rows),'compatible':sum(r['decision_compatible'] for r in rows),'counts':counts,'records_sha256':sha(a.out/'records.json'),'references_sha256':sha(a.out/'references.json'),'config_sha256':sha(a.out/'locked_config.json'),'failures':[{'network':r['network'],'degree':r['degree'],'mesh':r['intervals']} for r in rows if not r['reference_converged'] or not r['reference_overlap'] or not r['decision_compatible']]}
 if len(rows)!=summary['expected']:raise ValueError('incorrect number of records')
 atom(a.out/'summary.json',summary);print(json.dumps(summary),flush=True);return int(bool(summary['failures']))
if __name__=='__main__':sys.exit(main())
