"""Recheck a published frozen polynomial PINN fixture without prior run outputs.

Example: python scripts/106_portable_polynomial_check_20261101.py \
 --fixture configs/frozen_poly_new_training_2026-11-01.json --out /tmp/poly-check

Research prototype: interval assertions depend on the mathematical hypotheses
and numerical backends. The output directory must not exist. No PyTorch needed.
"""
from pathlib import Path
from fractions import Fraction
from collections import Counter
import argparse,hashlib,json,os,sys,traceback
from arb_certification import exact_fraction
from arb_cubic_reconstruction import CubicVerifier
from arb_shaped_remainder import ShapedVerifier
ROOT=Path(__file__).resolve().parents[1]
def write_json(path,value):
 data=(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n').encode()
 tmp=path.with_name(path.name+'.tmp')
 with tmp.open('xb') as f:f.write(data);f.flush();os.fsync(f.fileno())
 os.replace(tmp,path)
def run():
 parser=argparse.ArgumentParser();parser.add_argument('--fixture',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);parser.add_argument('--intervals',type=int,default=64);args=parser.parse_args()
 if args.intervals<1:parser.error('intervals must be positive')
 fixture=args.fixture.resolve();raw=fixture.read_bytes();data=json.loads(raw);nets=data['networks'];sha=hashlib.sha256(raw).hexdigest()
 if not nets or len({n['label'] for n in nets})!=len(nets):raise ValueError('fixture must have unique nonempty network labels')
 out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
 # Reuse the existing prototype definitions without executing its original driver.
 source=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0]
 original="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
 if source.count(original)!=1:raise ValueError('prototype driver layout changed')
 scope={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(source.replace(original,'OUT=None'),'polynomial_prototype_definitions','exec'),scope)
 config={'fixture_sha256':sha,'networks':len(nets),'intervals':args.intervals,'degrees':[2,3],'reference':'mpmath interval reference computed after both certificate bounds per network; rel 1e-7 abs 1e-14 cap 5000','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'prototype_sha256':hashlib.sha256((ROOT/'scripts/74_polynomial_forcing_prototype.py').read_bytes()).hexdigest(),'scope':'frozen binary64 weights as exact real constants, conditional on interval backend and assumptions; not public CLI'};write_json(out/'config.json',config)
 rows=[];refs={}
 for net in nets:
  name=net['label'];pair=[]
  for degree in (2,3):
   p=scope['PolynomialForcing'](fixture,name,sha,net['source_coefficients']);v=ShapedVerifier(p,'moment') if degree==2 else CubicVerifier(p,'moment');v.uniform(args.intervals);r=v.certify();l,u=exact_fraction(r['lower']),exact_fraction(r['upper']);assert 0<=l<=u
   row={'network':name,'degree':degree,'intervals':args.intervals,'lower':str(l),'upper':str(u),'remainder_upper':str(exact_fraction(r['remainder_upper']))};rows.append(row);pair.append(row)
  ref=scope['PolynomialReference'](net,net['source_coefficients']).maximum(relative=Fraction(1,10**7),absolute=Fraction(1,10**14),cap=5000);refs[name]=ref;rl,ru=scope['fr'](ref['lower']),scope['fr'](ref['upper'])
  for row in pair:
   l,u=Fraction(row['lower']),Fraction(row['upper']);row.update(reference_converged=ref['converged'],reference_overlap=l<=ru and u>=rl,reference_contained=l<=rl and u>=ru)
  write_json(out/'records.json',rows);write_json(out/'references.json',refs);print(name,'completed',flush=True)
 # Verify persisted records and exact set of requested network/degree keys.
 saved=json.loads((out/'records.json').read_text());expected={(n['label'],d) for n in nets for d in (2,3)};keys=[(r['network'],r['degree']) for r in saved]
 if len(keys)!=len(expected) or set(keys)!=expected:raise ValueError('persisted records incomplete or duplicate')
 b=(out/'records.json').read_bytes();summary={'records':len(saved),'expected':len(expected),'records_sha256':hashlib.sha256(b).hexdigest(),'reference_converged':sum(x['converged'] for x in refs.values()),'references':len(refs),'overlap':sum(x['reference_overlap'] for x in saved),'contained':sum(x['reference_contained'] for x in saved),'failures':[{'network':x['network'],'degree':x['degree']} for x in saved if not x['reference_overlap'] or not x['reference_converged']]};write_json(out/'summary.json',summary);print(json.dumps(summary),flush=True);return int(bool(summary['failures']))
if __name__=='__main__':sys.exit(run())
