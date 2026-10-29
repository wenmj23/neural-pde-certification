"""Off-center derivative diagnostics, sixth activation identity and input guards."""
from datetime import datetime,timezone
from fractions import Fraction
import hashlib
import importlib.util
import json
from math import factorial
from pathlib import Path
import mpmath as mp
from flint import arb,ctx
from arb_deep_network import DeepPoisson,activate
from arb_deep_refinements import RefinedDeepPoisson

ROOT=Path(__file__).resolve().parents[1]
FIXTURE=ROOT/'configs/frozen_networks_depth_2026-10-29.json'
DRIVER=ROOT/'scripts/23_depth_derivative_checks.py'
spec=importlib.util.spec_from_file_location('checks23',DRIVER)
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)


def main():
    out=ROOT/'results/generated/26_depth_interior_audit'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True,exist_ok=False);mp.mp.dps=90;ctx.prec=103
    rows=json.loads(FIXTURE.read_text())['networks'];sha=hashlib.sha256(FIXTURE.read_bytes()).hexdigest();checks=[]
    sources=[Path(__file__),DRIVER,ROOT/'scripts/arb_deep_network.py',ROOT/'scripts/arb_deep_refinements.py',FIXTURE]
    config=dict(status='running',scope='Off-center high-precision diagnostics and guards, not formal software audit',sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    def save(name,value):(out/name).write_text(json.dumps(value,indent=2)+'\n')
    save('config.json',config)
    def check(name,passed):
        checks.append(dict(name=name,passed=bool(passed)));save('checks.json',checks)
        if not passed:raise AssertionError(name)
    t5=[16,0,-136,0,240,0,-120];t6=[0]*8
    for j in range(1,len(t5)):
        t6[j-1]+=j*t5[j];t6[j+1]-=j*t5[j]
    check('sixth_polynomial_exact_recurrence',t6==[0,-272,0,1232,0,-1680,0,720])
    # Nonlinear inner jets exercise composition, not only tanh of an affine x.
    for i,coeff in enumerate([[Fraction(1,4),Fraction(2),Fraction(-1,3),Fraction(1,7),Fraction(2,9),Fraction(1,11),Fraction(-1,13)],
                              [Fraction(-3,5),Fraction(-1),Fraction(1,2),Fraction(0),Fraction(1,5),Fraction(0),Fraction(1,8)]]):
        mpc=[mp.mpf(q.numerator)/q.denominator for q in coeff]
        reference=mp.taylor(lambda x:mp.tanh(sum(c*x**j for j,c in enumerate(mpc))),mp.mpf(0),6)
        values=activate([arb(q.numerator)/q.denominator for q in coeff],'chain')
        for j,value in enumerate(values):check(f'sixth_composition_{i}_{j}',value.overlaps(arb(str(reference[j]))))
    z=[arb(-25).union(arb(88)),arb(1),arb(0),arb(0),arb(0),arb(0)]
    check('wide_monotone_fallback_finite',all(x.is_finite() for x in activate(z,'chain')))
    for row in rows:
        function=old.mp_function(row);omega=row['frequency']*mp.pi
        for a,b in [(Fraction(1,7),Fraction(2,7)),(Fraction(3,7),Fraction(4,7)),(Fraction(5,7),Fraction(6,7))]:
            points=[(3*a+b)/4,(a+3*b)/4];refs=[]
            for point in points:
                x=mp.mpf(point.numerator)/point.denominator;coeff=mp.taylor(function,x,5);d=[coeff[j]*factorial(j) for j in range(6)];q=x*(1-x);qp=1-2*x
                refs.append((-q*d[2]-2*qp*d[1]+2*d[0]-omega**2*mp.sin(omega*x),
                             -q*d[4]-4*qp*d[3]+12*d[2]+omega**4*mp.sin(omega*x),
                             -q*d[5]-5*qp*d[4]+20*d[3]+omega**5*mp.cos(omega*x)))
            cell=(arb(a.numerator)/a.denominator).lower().union((arb(b.numerator)/b.denominator).upper())
            for method in ('exp','chain','centered','value_centered','centered_lite','split'):
                cls=DeepPoisson if method in ('exp','chain','centered') else RefinedDeepPoisson
                p=cls(FIXTURE,row['label'],sha,method)
                enclosures=[p.residual(cell),p.r_second(cell),p.r_third(cell)]
                for i,values in enumerate(refs):
                    for j,(enclosure,value) in enumerate(zip(enclosures,values)):
                        check(f'interior_{row["label"]}_{method}_{a}_{i}_{j}',enclosure.overlaps(arb(str(value))))
        print('INTERIOR',row['label'],flush=True)
    # Invalid metadata is rejected even if its supplied checksum matches.
    fixture=json.loads(FIXTURE.read_text())
    for name,change in [('nonfinite_weight',lambda d:d['networks'][0]['layers'][0]['weight'][0].__setitem__(0,float('nan'))),
                        ('bad_shape',lambda d:d['networks'][0]['layers'][0]['weight'][0].append(1.)),
                        ('bad_activation',lambda d:d['networks'][0]['layers'][0].__setitem__('activation','relu'))]:
        modified=json.loads(json.dumps(fixture));change(modified);path=out/(name+'.json');path.write_text(json.dumps(modified));h=hashlib.sha256(path.read_bytes()).hexdigest();rejected=False
        try:DeepPoisson(path,rows[0]['label'],h)
        except ValueError:rejected=True
        check('reject_'+name,rejected)
    config.update(status='complete',checks_passed=len(checks));save('config.json',config)
    print('COMPLETED',out,'checks',len(checks),flush=True)


if __name__=='__main__':main()
