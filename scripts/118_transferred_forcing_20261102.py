"""Fixed frozen networks on three new polynomial forcings; no retraining.

f=12x²,20x³,42x⁵ with u=x-x⁴,x-x⁵,x-x⁷ respectively.
Compute degree-2/3 certificates first, then independent mpmath interval references.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,os,sys
from arb_certification import exact_fraction
from arb_cubic_reconstruction import CubicVerifier
from arb_shaped_remainder import ShapedVerifier
ROOT=Path(__file__).resolve().parents[1]
SOURCES={'quadratic_12x2':[0,0,12],'cubic_20x3':[0,0,0,20],'quintic_42x5':[0,0,0,0,0,42]}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atomic(p,obj):
 data=(json.dumps(obj,indent=2,allow_nan=False)+'\n').encode();tmp=p.with_suffix('.tmp')
 with tmp.open('xb') as f:f.write(data);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--fixtures',type=Path,nargs='+',required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 src=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
 if src.count(needle)!=1:raise ValueError('prototype layout changed')
 ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(src.replace(needle,'OUT=None'),'polydefs','exec'),ns)
 fixtures=[];rows=[];references=[]
 for fixture in a.fixtures:
  digest=sha(fixture);nets=json.loads(fixture.read_bytes())['networks'];fixtures.append({'sha256':digest,'networks':len(nets)})
  if len({n['label'] for n in nets})!=len(nets):raise ValueError('duplicate network')
  for net in nets:
   for name,coeff in SOURCES.items():
    pending=[]
    for degree in (2,3):
     p=ns['PolynomialForcing'](fixture,net['label'],digest,coeff);v=ShapedVerifier(p,'moment') if degree==2 else CubicVerifier(p,'moment');v.uniform(64);r=v.certify();lo,hi=exact_fraction(r['lower']),exact_fraction(r['upper'])
     if not 0<=lo<=hi:raise ValueError('invalid certificate')
     pending.append({'fixture_sha256':digest,'network':net['label'],'source':name,'degree':degree,'intervals':64,'lower':str(lo),'upper':str(hi),'remainder_upper':str(exact_fraction(r['remainder_upper'])),'gap_1_02':lo>0 and hi<=F(51,50)*lo})
    ref=ns['PolynomialReference'](net,coeff).maximum(relative=F(1,10**6),absolute=F(1,10**12),cap=2000)
    rl,ru=ns['fr'](ref['lower']),ns['fr'](ref['upper']);converged=bool(ref['converged'])
    references.append({'fixture_sha256':digest,'network':net['label'],'source':name,'lower':str(rl),'upper':str(ru),'converged':converged,'splits':ref['splits']})
    for row in pending:
     lo,hi=F(row['lower']),F(row['upper']);row.update(reference_converged=converged,reference_overlap=lo<=ru and rl<=hi,reference_contained=lo<=rl and ru<=hi);rows.append(row)
   print(net['label'],'all sources done',flush=True)
 atomic(a.out/'records.json',rows);atomic(a.out/'references.json',references)
 n=sum(f['networks'] for f in fixtures);summary={'fixtures':fixtures,'sources':SOURCES,'intervals':64,'degrees':[2,3],'expected_records':n*len(SOURCES)*2,'records':len(rows),'expected_references':n*len(SOURCES),'references':len(references),'reference_converged':sum(x['converged'] for x in references),'reference_overlap':sum(x['reference_overlap'] for x in rows),'reference_contained':sum(x['reference_contained'] for x in rows),'gap_1_02':{str(d):sum(x['gap_1_02'] for x in rows if x['degree']==d) for d in (2,3)},'records_sha256':sha(a.out/'records.json'),'references_sha256':sha(a.out/'references.json'),'script_sha256':sha(Path(__file__)),'failures':[{'network':x['network'],'source':x['source'],'degree':x['degree']} for x in rows if not x['reference_converged'] or not x['reference_overlap']],'scope':'transfer frozen weights to new polynomial forcing without retraining; conditional independent continuous reference; not public CLI'}
 if len(rows)!=summary['expected_records'] or len(references)!=summary['expected_references']:raise ValueError('missing rows')
 atomic(a.out/'summary.json',summary);print(json.dumps(summary),flush=True);return int(bool(summary['failures']))
if __name__=='__main__':sys.exit(main())
