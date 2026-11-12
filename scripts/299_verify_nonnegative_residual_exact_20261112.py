"""Independent symbolic post-exit check of nonnegative-residual controls."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'results/2026-11-12/nonnegative_residual_exact_controls'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rat(t):return s.Rational(t.numerator,t.denominator)
def main():
    summary=json.loads((D/'summary.json').read_text());rows=json.loads((D/'records.json').read_text())
    assert sha(D/'records.json')==summary['records_sha256']
    assert sha(ROOT/'scripts/298_nonnegative_residual_exact_controls_20261112.py')==summary['script_sha256']
    assert len(rows)==7==summary['expected']
    t=s.symbols('t',real=True)
    for row in rows:
        h=rat(F(row['support_width']));x=rat(F(row['maximum_location']));E=rat(F(row['true_max_error']))
        slope=h/2-h*h/6
        left=slope*t-t*t/2+t**3/(6*h);right=h*h*(1-t)/6
        assert s.simplify(-s.diff(left,t,2)-(1-t/h))==0
        assert s.diff(right,t,2)==0
        assert left.subs(t,0)==right.subs(t,1)==0
        assert s.simplify(left.subs(t,h)-right.subs(t,h))==0
        assert s.simplify(s.diff(left,t).subs(t,h)-s.diff(right,t).subs(t,h))==0
        assert s.diff(left,t).subs(t,x)==0 and left.subs(t,x)==E
        assert 0<x<h<1 and slope>0 and s.diff(right,t)<0
        assert F(row['global_over_true_max'])==F(1,8)/F(row['true_max_error'])
    out={'records':len(rows),'symbolic_equation_and_boundary_checks':len(rows),
         'ratio_first':rows[0]['global_over_true_max_approx'],
         'ratio_last':rows[-1]['global_over_true_max_approx'],
         'readback_script_sha256':sha(Path(__file__)),
         'scope':'exact symbolic algebra and rational readback for artificial piecewise-polynomial residuals; not a trained-network validation'}
    p=D/'independent_readback.json'
    if p.exists():raise FileExistsError(p)
    p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
