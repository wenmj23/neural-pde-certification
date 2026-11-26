"""Exact rational controls for a newly separated research forcing adapter."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import argparse, hashlib, json
from flint import arb
from arb_polynomial_forcing import PolynomialForcing, coefficient
from arb_reconstruction import point
from arb_certification import exact_fraction
from arb_shaped_remainder import ShapedVerifier
from arb_cubic_reconstruction import CubicVerifier

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    fixture=ROOT/'configs/frozen_poly_new_training_2026-11-01.json';h=sha(fixture);net=json.loads(fixture.read_text())['networks'][0]
    src=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0]
    needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
    assert src.count(needle)==1
    ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')}
    exec(compile(src.replace(needle,'OUT=None'),'legacy definitions','exec'),ns)
    coefficients=[[2],[0,0,12],[F(1,2)],[F(-1,2)],[F(1,2),F(-3,4),F(5,8)],[F(1,3),0,0,0,0,F(-7,11)]]
    rows=[]
    for coeff in coefficients:
        old=ns['PolynomialForcing'](fixture,net['label'],h,coeff)
        new=PolynomialForcing(fixture,net['label'],h,coeff)
        for j in range(5):
            for x in [F(0),F(1,7),F(1,2),F(6,7),F(1)]:
                target=sum((F(c)*F(factorial(i),factorial(i-j))*x**(i-j) for i,c in enumerate(coeff) if i>=j),F(0))
                ov,nv=old.source_derivative(point(x),j),new.source_derivative(point(x),j)
                ol,ou=exact_fraction(ov.lower()),exact_fraction(ov.upper());nl,nu=exact_fraction(nv.lower()),exact_fraction(nv.upper())
                assert nl<=target<=nu
                rows.append({'coefficients':[str(c) for c in coeff],'order':j,'point':str(x),'exact':str(target),'new_lower':str(nl),'new_upper':str(nu),'legacy_lower':str(ol),'legacy_upper':str(ou),'legacy_contains':ol<=target<=ou})
    rejected=0
    for val in [True,None,float('nan'),float('inf'),'1/0','nan',{}]:
        try:coefficient(val)
        except ValueError:rejected+=1
        else:raise AssertionError('invalid coefficient accepted')
    assert coefficient(0.1)==F.from_float(0.1) and coefficient('0.1')==F(1,10)
    comparisons=[]
    for raw in json.loads(fixture.read_text())['networks']:
        for deg in (2,3):
            pair=[]
            for cls in [ns['PolynomialForcing'],PolynomialForcing]:
                p=cls(fixture,raw['label'],h,raw['source_coefficients']);v=ShapedVerifier(p,'moment') if deg==2 else CubicVerifier(p,'moment');v.uniform(64);r=v.certify()
                pair.append({k:str(exact_fraction(r[k])) for k in ['lower','upper','remainder_upper']})
            assert pair[0]==pair[1]
            comparisons.append({'network':raw['label'],'degree':deg,'old':pair[0],'new':pair[1],'exact_match':True})
    (a.out/'source_controls.json').write_text(json.dumps(rows,indent=2)+'\n');(a.out/'integer_compatibility.json').write_text(json.dumps(comparisons,indent=2)+'\n')
    s={'rational_controls':len(rows),'new_contains':len(rows),'legacy_misses':sum(not r['legacy_contains'] for r in rows),'integer_certificates_equal':len(comparisons),'invalid_inputs_rejected':rejected,'fixture_sha256':h,'source_sha256':{p:sha(ROOT/p) for p in ['scripts/74_polynomial_forcing_prototype.py','scripts/arb_polynomial_forcing.py','scripts/316_audit_polynomial_coefficients_20261112.py']},'files':{p.name:sha(p) for p in a.out.glob('*.json')},'scope':'Legacy adapter floor-divides noninteger coefficients; published integer cases are unchanged in this comparison. New adapter is research-only; public CLI untouched.'}
    (a.out/'summary.json').write_text(json.dumps(s,indent=2)+'\n');print(json.dumps(s))
if __name__=='__main__':main()
