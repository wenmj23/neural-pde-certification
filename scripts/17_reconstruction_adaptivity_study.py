"""Higher-order, local refinement and reuse study; frozen or analytic inputs.

Compare equal residual-sample budgets and complete request costs. New quintic
network derivative and quadratic residual remainder are tested independently.
Known solutions are used for checks only, never refinement or stopping.
"""
import argparse
from datetime import datetime,timezone
from fractions import Fraction
import gc
import hashlib
import json
from pathlib import Path
import statistics
import sys
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpmath as mp
import flint
from flint import arb,ctx
from arb_certification import finite,exact_fraction,export,bounds
from arb_problem_families import FrozenSinePoisson,SineError,BumpError,QuarticError
from arb_reconstruction import Verifier,norm_enclosure,split_bernstein

ROOT=Path(__file__).resolve().parents[1]
FIXTURE=ROOT/'configs/frozen_networks_round07.json'
FIXTURE_SHA=hashlib.sha256(FIXTURE.read_bytes()).hexdigest()


def frac(e):return Fraction(int(e['numerator']),int(e['denominator']))


def serial(r):
    return dict(**bounds(r),intervals=r['intervals'],degree=r['degree'],subdivision=r['subdivision'],
        z_lower=export(r['z_lower'],'down'),z_upper=export(r['z_upper'],'up'),
        remainder_upper=export(r['remainder_upper'],'up'),residual_evaluations=r['residual_evaluations'],
        curvature_evaluations=r['curvature_evaluations'],polynomial_cell_visits=r['polynomial_cell_visits'],
        elapsed_s=r['elapsed_s'])


def factories():
    fs={f'sine_{k}':(lambda k=k:SineError(k)) for k in (1,20,50)}
    fs['mixed_weak_high']=lambda:SineError(1,frequencies=[1,5,81],
        coefficients=[Fraction(120,128),Fraction(7,128),Fraction(1,128)],name='mixed_weak_high')
    fs['bump_0.02']=lambda:BumpError('0.02')
    fs['bump_0.005']=lambda:BumpError('0.005')
    fs['quartic']=QuarticError
    for k in (1,5):
        for seed in (7,17,27):
            label=f'PINN_k{k}_seed{seed}'
            fs[label]=lambda label=label:FrozenSinePoisson(FIXTURE,label,FIXTURE_SHA)
    return fs


def network_reference(problem,n=4096):
    """Different continuous formula at 169 bits, sharing the Arb backend."""
    nodes=[Fraction(i,n) for i in range(n+1)]
    values=[finite(problem.error(arb(i)/n),'reference error') for i in range(n+1)]
    lo=max(v.abs_lower() for v in values);hi=arb(0)
    for i,(a,b) in enumerate(zip(nodes[:-1],nodes[1:])):
        aa,bb=arb(a.numerator)/a.denominator,arb(b.numerator)/b.denominator
        K=finite(problem.residual(aa.lower().union(bb.upper())),'reference curvature').abs_upper();h=bb-aa
        hi=max(hi,(max(values[i].abs_upper(),values[i+1].abs_upper())+K*h*h/8).upper())
    return dict(lower=export(lo,'down'),upper=export(hi,'up'),intervals=n,bits=problem.bits,
                scope='continuous reference with known sin(k*pi*x); different formula, same Arb backend')


def request(factory,degree,adaptive=False,subdivision=2,target=Fraction(6,5),cap=16384,cache=True):
    gc.collect();start=time.perf_counter();cpu=time.process_time();p=factory()
    v=Verifier(p,degree,subdivision,cache=cache);n=32;history=[];success=False
    while n<=cap:
        if adaptive:v.refine_to(n)
        else:v.uniform(n)
        r=v.certify();s=serial(r);low,up=frac(s['lower']),frac(s['upper'])
        success=low>0 and up<=target*low
        history.append(dict(intervals=n,lower=s['lower'],upper=s['upper']))
        if success or n==cap:break
        n*=2
    return dict(**s,adaptive=adaptive,cache=cache,target_numerator=target.numerator,
        target_denominator=target.denominator,success=bool(success),history=history,
        wall_ms=(time.perf_counter()-start)*1000,cpu_ms=(time.process_time()-cpu)*1000)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path);parser.add_argument('--repeats',type=int,default=3)
    args=parser.parse_args()
    if args.repeats<3:raise ValueError('At least three timing repeats required')
    out=args.out or ROOT/'results/generated/17_reconstruction_adaptivity'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True,exist_ok=False);fs=factories();checks=[];references={}
    def check(name,ok):
        checks.append(dict(name=name,passed=bool(ok)))
        (out/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
        if not ok:raise AssertionError(name)
    # Exact T5 polynomial recurrence, independent of numerical evaluation.
    t4=[0,16,0,-40,0,24]
    derivative=[(i+1)*t4[i+1] for i in range(len(t4)-1)];t5=[0]*7
    for i,c in enumerate(derivative):t5[i]+=c;t5[i+2]-=c
    check('exact_fifth_tanh_derivative',t5==[16,0,-136,0,240,0,-120])
    # Product extremum for nodes 0,1/2,1: max |t(t-1/2)(t-1)|=1/(12sqrt3).
    ctx.prec=103;sqrt3=arb(3).sqrt();maximum=arb(1)/(12*sqrt3)
    for root in [(3-sqrt3)/6,(3+sqrt3)/6]:
        value=(root*(root-arb(1)/2)*(root-1)).abs_upper()
        check('quadratic_remainder_constant_'+str(len(checks)),value>=maximum.lower()
              and exact_fraction(value)-exact_fraction(maximum.upper())<Fraction(1,10**25))
    # Exact quadratic residual and subdivision: this isolates integration/extrema.
    for degree,depth in [(1,0),(2,0),(2,1),(2,2)]:
        p=QuarticError();v=Verifier(p,degree,depth);v.uniform(1);r=v.certify()
        check(f'quartic_known_norm_{degree}_{depth}',exact_fraction(r['lower'])<=p.known_error<=exact_fraction(r['upper']))
        if degree==2:
            check('zero_third_derivative_zero_remainder_'+str(depth),exact_fraction(r['remainder_upper'])==0)
            if depth:
                check('subdivision_recovers_quartic_peak_'+str(depth),
                      exact_fraction(r['upper'])-exact_fraction(r['lower'])<Fraction(1,10**25))
    # Arbitrary nonuniform rational mesh.
    p=QuarticError();v=Verifier(p,2,2)
    r=v.certify([Fraction(0),Fraction(1,7),Fraction(3,10),Fraction(4,5),Fraction(1)])
    check('nonuniform_quartic_norm',exact_fraction(r['lower'])<=p.known_error<=exact_fraction(r['upper']))
    for degree in (1,2):
        p=SineError(200);v=Verifier(p,degree,1);v.uniform(100);r=v.certify()
        check('aliased_samples_do_not_certify_zero_'+str(degree),exact_fraction(r['upper'])>=p.known_error
              and exact_fraction(r['remainder_upper'])>0)
    class Bad:
        boundary_error=(0,0);bits=103
        def guard(self):ctx.prec=103
        def residual(self,x):return arb(0)
        def r_second(self,x):return arb('nan')
        def r_third(self,x):return arb('nan')
    class Zero(Bad):
        def r_second(self,x):return arb(0)
        def r_third(self,x):return arb(0)
    for degree in (1,2):
        zero=Verifier(Zero(),degree,2).certify()
        check('zero_error_exact_bounds_'+str(degree),exact_fraction(zero['lower'])==exact_fraction(zero['upper'])==0)
    for degree in (1,2):
        try:Verifier(Bad(),degree).certify()
        except ArithmeticError:rejected=True
        else:rejected=False
        check('reject_nonfinite_curvature_'+str(degree),rejected)
    # Reusing nodes/cells must not alter the final interval on a nested mesh.
    comparison=[]
    for cache in (False,True):
        v=Verifier(SineError(20),2,1,cache=cache)
        for n in (16,32,64):v.uniform(n);r=v.certify()
        comparison.append(serial(r))
    check('cache_preserves_exact_final_bounds',all(frac(comparison[0][k])==frac(comparison[1][k]) for k in ['lower','upper']))
    check('cache_reduces_residual_evaluations',comparison[1]['residual_evaluations']<comparison[0]['residual_evaluations'])
    (out/'cache_sanity.json').write_text(json.dumps(comparison,indent=2)+'\n')
    print('Interpolation/integration, nonuniform mesh, aliasing, guards and cache tests passed.',flush=True)

    # Fifth-order derivative point/whole-cell diagnostics for every family.
    with mp.workdps(100):
        for name,factory in fs.items():
            p=factory();error=p.mp_function()
            for numerator,denominator in [(1,10),(3,10),(1,2),(7,10),(9,10)]:
                x=mp.mpf(numerator)/denominator;xx=arb(numerator)/denominator
                for order,method in [(2,p.residual),(4,p.r_second),(5,p.r_third)]:
                    value=-mp.diff(error,x,order)
                    # Numerical differences can leave a tiny nonzero remainder
                    # at exact odd-symmetry zeros; these zeros are analytic.
                    if order==5 and (isinstance(p,QuarticError) or (isinstance(p,BumpError) and 2*numerator==denominator)):
                        value=mp.mpf(0)
                    q=Fraction(*mp.libmp.to_rational(value._mpf_));enclosure=finite(method(xx),'derivative diagnostic')
                    check(f'derivative_point_{name}_{order}_{numerator}',exact_fraction(enclosure.lower())<=q<=exact_fraction(enclosure.upper()))
            for i in (0,3,7):
                cell=(arb(i)/8).union(arb(i+1)/8);enclosure=finite(p.r_third(cell),'third derivative cell')
                for j in (0,1,2):
                    x=mp.mpf(2*i+j)/16;value=-mp.diff(error,x,5)
                    if isinstance(p,QuarticError) or (isinstance(p,BumpError) and x==mp.mpf(1)/2):value=mp.mpf(0)
                    q=Fraction(*mp.libmp.to_rational(value._mpf_))
                    check(f'derivative_cell_{name}_{i}_{j}',exact_fraction(enclosure.lower())<=q<=exact_fraction(enclosure.upper()))
            print('Derivative checks passed:',name,flush=True)
    for name,factory in fs.items():
        p=factory()
        if hasattr(p,'known_error'):
            references[name]=dict(kind='analytic_exact_maximum',numerator=str(p.known_error.numerator),denominator=str(p.known_error.denominator))
        else:
            p=FrozenSinePoisson(FIXTURE,name,FIXTURE_SHA,169)
            references[name]=network_reference(p)
            print('Continuous reference computed:',name,flush=True)
    (out/'references.json').write_text(json.dumps(references,indent=2)+'\n')
    def consistent(name,r,label):
        lo,hi=frac(r['lower']),frac(r['upper']);ref=references[name]
        if ref.get('kind')=='analytic_exact_maximum':
            exact=Fraction(int(ref['numerator']),int(ref['denominator']));ok=lo<=exact<=hi
        else:
            # A newly tight certificate can be narrower than the reference.
            # Require compatibility, not containment of a wider reference.
            ok=lo<=frac(ref['upper']) and hi>=frac(ref['lower'])
        check(label,ok)
    screen=[]
    for name,factory in fs.items():
        for budget in (129,513,1025):
            for degree,depth in [(1,0),(2,0),(2,2)]:
                p=factory();v=Verifier(p,degree,depth);n=(budget-1)//degree;v.uniform(n)
                r=serial(v.certify());r.update(problem=name,sample_budget=budget)
                check(f'equal_budget_count_{name}_{budget}_{degree}_{depth}',r['residual_evaluations']==budget)
                consistent(name,r,f'screen_reference_{name}_{budget}_{degree}_{depth}')
                screen.append(r)
        (out/'equal_sample_screen.json').write_text(json.dumps(screen,indent=2)+'\n')
        print('Equal residual-sample budgets checked:',name,flush=True)
    # Finite, bounded matched-gap comparison. Report cap exhaustion explicitly.
    selected=['sine_20','sine_50','bump_0.005','PINN_k1_seed7','PINN_k5_seed7','PINN_k5_seed17']
    methods=[(1,False,2),(1,True,2),(2,False,2),(2,True,2)]
    runs=[]
    for repeat in range(args.repeats):
        for name in selected:
            ordered=methods[repeat%len(methods):]+methods[:repeat%len(methods)]
            for degree,adaptive,depth in ordered:
                r=request(fs[name],degree,adaptive,depth);r.update(problem=name,repeat=repeat)
                consistent(name,r,f'timed_reference_{name}_{repeat}_{degree}_{adaptive}')
                if r['success']:
                    check(f'timed_stop_{name}_{repeat}_{degree}_{adaptive}',frac(r['lower'])>0
                          and frac(r['upper'])<=Fraction(6,5)*frac(r['lower']))
                runs.append(r);(out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n')
                print(f"repeat={repeat+1} {name} degree={degree} adaptive={adaptive}: success={r['success']} n={r['intervals']} ms={r['wall_ms']:.2f}",flush=True)
    extra=[]
    for repeat in range(args.repeats):
        for degree,adaptive in [(1,False),(2,False),(2,True)]:
            r=request(fs['PINN_k1_seed7'],degree,adaptive,2,Fraction(501,500));r.update(problem='PINN_k1_seed7',repeat=repeat)
            consistent('PINN_k1_seed7',r,f'strict_main_reference_{repeat}_{degree}_{adaptive}')
            extra.append(r)
        for cache in (False,True):
            r=request(fs['bump_0.005'],2,True,2,cache=cache);r.update(problem='bump_0.005',repeat=repeat)
            consistent('bump_0.005',r,f'cache_cost_reference_{repeat}_{cache}');extra.append(r)
    (out/'extra_runs.json').write_text(json.dumps(extra,indent=2)+'\n')
    summary=[]
    for name in selected:
        for degree,adaptive,depth in methods:
            rows=[r for r in runs if r['problem']==name and r['degree']==degree and r['adaptive']==adaptive]
            check(f'deterministic_mesh_{name}_{degree}_{adaptive}',len({(r['intervals'],r['success']) for r in rows})==1)
            summary.append(dict(problem=name,degree=degree,adaptive=adaptive,subdivision=depth,
                success=rows[0]['success'],intervals=rows[0]['intervals'],residual_evaluations=rows[0]['residual_evaluations'],
                curvature_evaluations=rows[0]['curvature_evaluations'],
                median_wall_ms=statistics.median(r['wall_ms'] for r in rows),min_wall_ms=min(r['wall_ms'] for r in rows),
                max_wall_ms=max(r['wall_ms'] for r in rows)))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    # Absolute tolerance decisions; upper>tau alone never implies rejection.
    decisions=[]
    for name in ['PINN_k1_seed7','PINN_k5_seed7','PINN_k5_seed17']:
        row=next(r for r in runs if r['problem']==name and r['degree']==2 and not r['adaptive'])
        tau=Fraction(1,1000);lo,hi=frac(row['lower']),frac(row['upper'])
        status='accept' if hi<=tau else ('reject' if lo>tau else 'uncertain')
        ref=references[name]
        if status=='accept':check('reference_supports_accept_'+name,frac(ref['upper'])<=tau)
        elif status=='reject':check('reference_supports_reject_'+name,frac(ref['lower'])>tau)
        decisions.append(dict(problem=name,tolerance='1/1000',status=status,lower=row['lower'],upper=row['upper']))
    (out/'decisions.json').write_text(json.dumps(decisions,indent=2)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    for ax,name in zip(axes,['sine_50','bump_0.005']):
        for degree,depth,label in [(1,0,'Linear'),(2,0,'Quadratic'),(2,2,'Quadratic + norm subdivision')]:
            rows=[r for r in screen if r['problem']==name and r['degree']==degree and r['subdivision']==depth]
            ax.loglog([r['sample_budget'] for r in rows],[float(frac(r['upper'])/Fraction(1,10000)) for r in rows],'-o',label=label)
        ax.set(title=name,xlabel='Residual evaluation budget',ylabel='Upper bound / exact error');ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.tight_layout();fig.savefig(out/'equal_budget.png',dpi=160);plt.close(fig)
    sources=[ROOT/'scripts'/n for n in ['17_reconstruction_adaptivity_study.py','arb_problem_families.py','arb_reconstruction.py','arb_certification.py']]
    metadata=dict(date='2026-10-29',fixture_sha256=FIXTURE_SHA,networks='six frozen outputs, no retraining',
        manufactured='seven analytic error families with known continuous norms',bits=103,reference_bits=169,
        repeats=args.repeats,screen_configurations=len(screen),timed_requests=len(runs)+len(extra),checks_passed=len(checks),
        sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        python=sys.version,python_flint=flint.__version__,flint=flint.__FLINT_VERSION__,
        assumptions=['zero error boundary','linear residual requires r C2; quadratic requires r C3 (e.g. v C5, f C3)',
                     'Arb interval-function correctness','fixed oracle and precision for cache reuse'],
        timing='Fresh model and cache; all refinement attempts/indicators/integrals counted; gc.collect before timer. Research references and file writes excluded.',
        limitations=['Network reference uses a different formula at higher precision but shares Arb; checks require overlap if reference is wider',
                     'High-precision differentiation is diagnostic, not a formal proof','Assistant review only; no independent software audit',
                     'No claim that every higher-order or adaptive strategy is faster or tighter'])
    (out/'config.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print('FINAL checks',len(checks),'screen configs',len(screen),'timed requests',len(runs)+len(extra),flush=True)
    print('RESULT_DIRECTORY',out,flush=True)


if __name__=='__main__':main()
