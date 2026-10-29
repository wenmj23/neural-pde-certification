"""Cross-check derivative jets against high-precision mpmath differentiation.

Numerical differentiation is diagnostic, not a formal proof. Interval enclosure
construction relies on the analytic chain and mean-value identities in the code.
"""
from datetime import datetime,timezone
from fractions import Fraction
import hashlib
import json
from math import factorial
from pathlib import Path
import mpmath as mp
from flint import arb,ctx
from arb_deep_network import DeepPoisson
from arb_certification import exact_fraction,export

ROOT=Path(__file__).resolve().parents[1]
FIXTURE=ROOT/'configs/frozen_networks_depth_2026-10-29.json'


def mp_function(row):
    weights=[([[mp.mpf(float(w)) for w in r] for r in l['weight']],[mp.mpf(float(b)) for b in l['bias']]) for l in row['layers']]
    def function(x):
        state=[x]
        for i,(w,b) in enumerate(weights):
            state=[bb+sum(ww*v for ww,v in zip(ri,state)) for ri,bb in zip(w,b)]
            if i<len(weights)-1:state=[mp.tanh(v) for v in state]
        return state[0]
    return function


def main():
    out=ROOT/'results/generated/23_depth_derivatives'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True,exist_ok=False);mp.mp.dps=90
    sha=hashlib.sha256(FIXTURE.read_bytes()).hexdigest();rows=json.loads(FIXTURE.read_text())['networks'];checks=[];records=[]
    sources=[Path(__file__),ROOT/'scripts/arb_deep_network.py',ROOT/'scripts/arb_certification.py',ROOT/'scripts/arb_problem_families.py',FIXTURE]
    (out/'source_snapshot').mkdir()
    for p in sources:
        if p.suffix=='.py':(out/'source_snapshot'/p.name).write_bytes(p.read_bytes())
    config=dict(status='running',sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},precision_bits=103,reference_decimal_digits=90,scope='High-precision numerical differentiation and interval diagnostics, not formal verification')
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    def check(name,passed):
        checks.append(dict(name=name,passed=bool(passed)))
        (out/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
        if not passed:raise AssertionError(name)
    points=[Fraction(1,7),Fraction(1,2),Fraction(6,7)]
    for row in rows:
        function=mp_function(row);refs={}
        for x in points:refs[x]=mp.taylor(function,mp.mpf(x.numerator)/x.denominator,5)
        for method in ('exp','chain','centered'):
            p=DeepPoisson(FIXTURE,row['label'],sha,method)
            for x in points:
                jet=p.network_jet(arb(x.numerator)/x.denominator,5)
                for j,v in enumerate(jet):
                    reference=arb(str(refs[x][j]))
                    check(f'point_{row["label"]}_{method}_{x}_{j}',v.overlaps(reference))
            for left,right in [(Fraction(0),Fraction(1,4)),(Fraction(1,4),Fraction(3,4)),(Fraction(3,4),Fraction(1))]:
                cell=(arb(left.numerator)/left.denominator).lower().union((arb(right.numerator)/right.denominator).upper())
                try:jet=p.network_jet(cell,5)
                except ArithmeticError as error:
                    # Wide native exp series can be indeterminate; safe rejection
                    # is recorded, never converted into a zero bound.
                    records.append(dict(label=row['label'],method=method,cell=[str(left),str(right)],status='safe_rejection',reason=str(error)))
                    check(f'wide_rejection_only_exp_{row["label"]}_{method}_{left}',method=='exp')
                    continue
                x=(left+right)/2;ref=mp.taylor(function,mp.mpf(x.numerator)/x.denominator,5)
                for j,v in enumerate(jet):check(f'cell_{row["label"]}_{method}_{left}_{j}',v.overlaps(arb(str(ref[j]))))
                records.append(dict(label=row['label'],method=method,cell=[str(left),str(right)],status='finite',fifth_derivative_width=str((jet[5].upper()-jet[5].lower())*factorial(5))))
        print('CHECKED',row['label'],flush=True)
        (out/'records.json').write_text(json.dumps(records,indent=2)+'\n')
    p=DeepPoisson(FIXTURE,rows[0]['label'],sha);old_cap=ctx.cap;ctx.cap=old_cap+1;rejected=False
    try:p.residual(arb(1)/2)
    except ValueError:rejected=True
    ctx.cap=old_cap;check('series_cap_change_rejected',rejected)
    config.update(status='complete',checks_passed=len(checks),wide_safe_rejections=sum(r['status']=='safe_rejection' for r in records))
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    print('COMPLETED',out,'checks',len(checks),'safe_rejections',config['wide_safe_rejections'],flush=True)


if __name__=='__main__':main()
