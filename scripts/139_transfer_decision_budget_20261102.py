"""Compare fixed-mesh decision certificates for transferred source problems.

Fixed tolerances are post-hoc relative to previously computed references.
No wall-clock comparisons or claims of blind evaluation.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,os,sys
from arb_certification import exact_fraction
from arb_cubic_reconstruction import CubicVerifier
from arb_shaped_remainder import ShapedVerifier
ROOT=Path(__file__).resolve().parents[1]
SOURCES={'quadratic_12x2':[0,0,12],'cubic_20x3':[0,0,0,20],'quintic_42x5':[0,0,0,0,0,42]};MESH=(8,16,32,64);DEGREES=(2,3);TOL=(F(1,10),F(1,5),F(3,10))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2)+'\n').encode();tmp=p.with_suffix('.tmp')
 with tmp.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def classify(lo,hi,t):return 'within' if hi<=t else ('exceeds' if lo>t else 'undetermined')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--fixtures',nargs='+',required=True,type=Path);ap.add_argument('--references',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 prev=json.loads((a.references/'summary.json').read_text());rp=a.references/'references.json'
 if sha(rp)!=prev['references_sha256'] or prev['sources']!=SOURCES:raise ValueError('reference inputs mismatch')
 refs={(r['fixture_sha256'],r['network'],r['source']):r for r in json.loads(rp.read_text())};expected_refs={(sha(f),n['label'],source) for f in a.fixtures for n in json.loads(f.read_bytes())['networks'] for source in SOURCES}
 if set(refs)!=expected_refs or len(refs)!=len(expected_refs) or any(not r['converged'] for r in refs.values()):raise ValueError('reference keys or convergence mismatch')
 src=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
 if src.count(needle)!=1:raise ValueError('prototype changed')
 ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(src.replace(needle,'OUT=None'),'defs','exec'),ns)
 rows=[]
 for fixture in a.fixtures:
  h=sha(fixture);nets=json.loads(fixture.read_bytes())['networks']
  for net in nets:
   for name,coeff in SOURCES.items():
    ref=refs[(h,net['label'],name)];rl,ru=F(ref['lower']),F(ref['upper'])
    for degree in DEGREES:
     p=ns['PolynomialForcing'](fixture,net['label'],h,coeff);v=ShapedVerifier(p,'moment') if degree==2 else CubicVerifier(p,'moment')
     for mesh in MESH:
      v.uniform(mesh);r=v.certify();lo,hi=exact_fraction(r['lower']),exact_fraction(r['upper']);decisions={}
      for tol in TOL:
       dec=classify(lo,hi,tol);refdec=classify(rl,ru,tol);decisions[str(tol)]={'decision':dec,'reference_decision':refdec,'compatible':dec=='undetermined' or refdec=='undetermined' or dec==refdec}
      rows.append({'fixture_sha256':h,'network':net['label'],'source':name,'degree':degree,'intervals':mesh,'lower':str(lo),'upper':str(hi),'reference_lower':str(rl),'reference_upper':str(ru),'reference_overlap':lo<=ru and rl<=hi,'gap_1_02':lo>0 and hi<=F(51,50)*lo,'decisions':decisions})
   print(net['label'],'done',flush=True)
 atom(a.out/'records.json',rows);summary={'fixtures':prev['fixtures'],'sources':SOURCES,'meshes':MESH,'degrees':DEGREES,'tolerances':[str(t) for t in TOL],'expected':len(expected_refs)*len(MESH)*len(DEGREES),'records':len(rows),'reference_sha256':prev['references_sha256'],'records_sha256':sha(a.out/'records.json'),'script_sha256':sha(Path(__file__)),'counts':{},'failures':[]}
 for degree in DEGREES:
  for mesh in MESH:
   group=[r for r in rows if r['degree']==degree and r['intervals']==mesh]
   for tol in TOL:
    vals=[r['decisions'][str(tol)] for r in group];summary['counts'][f'{degree}/{mesh}/{tol}']={'within':sum(v['decision']=='within' for v in vals),'exceeds':sum(v['decision']=='exceeds' for v in vals),'undetermined':sum(v['decision']=='undetermined' for v in vals),'reference_undetermined':sum(v['reference_decision']=='undetermined' for v in vals),'incompatible':sum(not v['compatible'] for v in vals)}
 summary['failures']=[{'network':r['network'],'source':r['source'],'degree':r['degree'],'mesh':r['intervals']} for r in rows if not r['reference_overlap'] or any(not d['compatible'] for d in r['decisions'].values())]
 if len(rows)!=summary['expected']:raise ValueError('count mismatch')
 atom(a.out/'summary.json',summary);print(json.dumps(summary),flush=True);return int(bool(summary['failures']))
if __name__=='__main__':sys.exit(main())
