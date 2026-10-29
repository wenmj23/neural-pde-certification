"""Validate and benchmark shaped quadratic interpolation remainders.

Known solutions and continuous references are validation-only. Fixed-grid
comparisons isolate the envelope; complete fresh requests measure costs.
"""
from datetime import datetime,timezone
from fractions import Fraction
import gc
import hashlib
import importlib.util
import json
from pathlib import Path
import statistics
import time

import flint
from flint import arb,ctx
from arb_certification import finite,exact_fraction,export
from arb_problem_families import SmoothProblem,QuarticError
from arb_reconstruction import Verifier
from arb_shaped_remainder import ShapedVerifier,shaped_norm,KAPPA,power_to_bernstein,split_exact

ROOT=Path(__file__).resolve().parents[1]
DRIVER=ROOT/'scripts/17_reconstruction_adaptivity_study.py'
spec=importlib.util.spec_from_file_location('study17',DRIVER)
study=importlib.util.module_from_spec(spec);spec.loader.exec_module(study)
REFERENCE=ROOT/'results/generated/17_reconstruction_adaptivity/20261029T034336591977Z/independent_reference_20261029T034918489790Z/references.json'


def frac(e):return Fraction(int(e['numerator']),int(e['denominator']))
def as_arb(q):return arb(q.numerator)/q.denominator


def exact_green(cells,x,shaped=False):
    """Independent rational integration; shaped=True returns D(x)/sqrt(3)."""
    total=Fraction(0)
    for a,b,delta in cells:
        h=b-a;cuts={Fraction(0),Fraction(1)}
        if shaped:cuts.add(Fraction(1,2))
        if a<x<b:cuts.add((x-a)/h)
        cuts=sorted(cuts)
        for left,right in zip(cuts[:-1],cuts[1:]):
            midpoint=(left+right)/2
            if a+h*midpoint<=x:alpha,beta=(1-x)*a,(1-x)*h
            else:alpha,beta=x*(1-a),-x*h
            if shaped:
                sign=1 if midpoint<Fraction(1,2) else -1
                source=[Fraction(0),Fraction(6*sign),Fraction(-18*sign),Fraction(12*sign)]
            else:source=[Fraction(1)]
            integrand=[Fraction(0)]*(len(source)+1)
            for j,c in enumerate(source):integrand[j]+=alpha*c;integrand[j+1]+=beta*c
            total+=delta*h*sum(c*(right**(j+1)-left**(j+1))/(j+1) for j,c in enumerate(integrand))
    return total


def correction_value(t):
    t=min(t,1-t)
    return Fraction(3,16)*t*t-t**3+Fraction(3,2)*t**4-Fraction(3,5)*t**5


def serialize(r):
    s=study.serial(r)
    s.update(mode=r['mode'],indicator=r['indicator'],indicator_peak=r['indicator_peak'],
             constant_remainder_upper=export(r['constant_remainder_upper'],'up'),
             raw_remainder_upper=export(r['raw_remainder_upper'],'up'))
    return s


def request(factory,mode,adaptive,indicator,target,cap=16384):
    gc.collect();start=time.perf_counter();cpu=time.process_time()
    p=factory();v=ShapedVerifier(p,mode,indicator);n=32;history=[]
    while True:
        if adaptive:v.refine_to(n)
        else:v.uniform(n)
        s=serialize(v.certify());lo,hi=frac(s['lower']),frac(s['upper'])
        success=lo>0 and hi<=target*lo
        history.append(dict(intervals=n,lower=s['lower'],upper=s['upper']))
        if success or n==cap:break
        n=min(cap,2*n)
    return dict(**s,adaptive=adaptive,target=str(target),success=bool(success),history=history,
                wall_ms=(time.perf_counter()-start)*1000,cpu_ms=(time.process_time()-cpu)*1000)


class ZeroProblem(SmoothProblem):
    def __init__(self):self.setup(103)
    def residual(self,x):self.guard();return arb(0)
    def r_third(self,x):self.guard();return arb(0)


def main():
    out=ROOT/'results/generated/20_shaped_remainder'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True,exist_ok=False)
    source_files=[Path(__file__),ROOT/'scripts/arb_shaped_remainder.py',DRIVER,
                  ROOT/'scripts/arb_reconstruction.py',ROOT/'scripts/arb_problem_families.py',ROOT/'scripts/arb_certification.py']
    snapshots=out/'source_snapshot';snapshots.mkdir()
    for p in source_files:(snapshots/p.name).write_bytes(p.read_bytes())
    config=dict(status='running',date='2026-10-29',bits=103,repeats=3,
                backend=dict(python_flint=flint.__version__,flint=flint.__FLINT_VERSION__),
                kappa=str(KAPPA),kappa_proof='Maximum exact rational Bernstein control coefficient of J(t)/sqrt(3) after four midpoint subdivisions on each symmetric half.',
                sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files+[REFERENCE,ROOT/'configs/frozen_networks_round07.json']},
                assumptions=['Unit interval zero error boundary','Residual C3 with included third derivative bound','Interval backend correctness; fixed network and precision','Assistant review and diagnostic checks; no independent/formal audit'],
                timing='Fresh model/cache; gc.collect before clock; all failed budgets, adaptive scores, integration, norm bounds and exports included; reference checks and file writes excluded.')
    def save(name,value):(out/name).write_text(json.dumps(value,indent=2)+'\n')
    save('config.json',config);checks=[]
    def check(name,ok):
        checks.append(dict(name=name,passed=bool(ok)));save('checks.json',checks)
        if not ok:raise AssertionError(name)
    fs=study.factories();refs=json.loads(REFERENCE.read_text())
    def consistent(name,s,label):
        pknown=name in ('sine_1','sine_20','sine_50','mixed_weak_high','bump_0.02','bump_0.005','quartic')
        if pknown:lo=hi=Fraction(1,10000)
        else:lo,hi=frac(refs[name]['lower']),frac(refs[name]['upper'])
        check(label,frac(s['lower'])<=lo and frac(s['upper'])>=hi)
    ctx.prec=103
    check('kappa_positive',KAPPA>0)
    # Exact polynomial proof data: mass=3/8 (after removing sqrt(3)),
    # symmetry centroid=1/2, and the local correction identity at rational x.
    check('single_cell_shaped_green_midpoint',exact_green([(Fraction(0),Fraction(1),Fraction(1))],Fraction(1,2),True)==Fraction(7,160))
    synthetic=[(Fraction(0),Fraction(1,7),Fraction(3,2)),
               (Fraction(1,7),Fraction(3,5),Fraction(1,100)),
               (Fraction(3,5),Fraction(1),Fraction(7,3))]
    for i in range(65):
        x=Fraction(i,64);constant=exact_green(synthetic,x);shaped=exact_green(synthetic,x,True)
        local=Fraction(0)
        for a,b,d in synthetic:
            if a<x<b:local=d*(b-a)**2*correction_value((x-a)/(b-a))
        check('exact_moment_identity_'+str(i),shaped==Fraction(3,8)*constant+local)
        check('correction_bound_sample_'+str(i),correction_value(x)<=KAPPA)
    # Check every rational Bernstein coefficient used in the correction proof.
    poly=[Fraction(0),Fraction(0),Fraction(3,64),Fraction(-1,8),Fraction(3,32),Fraction(-3,160)]
    pieces=[power_to_bernstein(poly)]
    for _ in range(4):pieces=[child for p in pieces for child in split_exact(p)]
    check('all_correction_control_coefficients_bounded',all(c<=KAPPA for p in pieces for c in p))
    exact_cell=[dict(a=arb(a.numerator)/a.denominator,b=as_arb(b),h=as_arb(b-a),delta=as_arb(d)) for a,b,d in synthetic]
    low,high=shaped_norm(exact_cell)
    for i in range(33):
        exact=as_arb(exact_green(synthetic,Fraction(i,32),True))*arb(3).sqrt()
        check('quintic_norm_covers_independent_point_'+str(i),exact_fraction(high)>=exact_fraction(exact.lower()))
    midpoint=7*arb(3).sqrt()/160
    lo,hi=shaped_norm([dict(a=arb(0),b=arb(1),h=arb(1),delta=arb(1))])
    check('symmetric_shaped_maximum',lo<=midpoint.upper() and hi>=midpoint.lower())
    for factory in [ZeroProblem,QuarticError]:
        for mode in ('constant','full','moment'):
            v=ShapedVerifier(factory(),mode);v.uniform(16);s=serialize(v.certify())
            check('zero_remainder_'+factory.__name__+'_'+mode,frac(s['remainder_upper'])==0)
    class InvalidProblem(ZeroProblem):
        def r_third(self,x):return arb('nan')
    for mode in ('constant','full','moment'):
        rejected=False
        try:ShapedVerifier(InvalidProblem(),mode).certify()
        except ArithmeticError:rejected=True
        check('reject_nan_'+mode,rejected)
    v=ShapedVerifier(ZeroProblem());ctx.prec=104;rejected=False
    try:v.certify()
    except ValueError:rejected=True
    ctx.prec=103;check('precision_guard',rejected)
    # Fixed meshes: same z and same calls; only D is changed.
    screen=[]
    for name,factory in fs.items():
        for budget in (129,513,1025):
            previous=None
            for mode in ('constant','full','moment'):
                p=factory();v=ShapedVerifier(p,mode);v.uniform((budget-1)//2);r=v.certify();s=serialize(r)
                s.update(problem=name,sample_budget=budget)
                check(f'sample_count_{name}_{budget}_{mode}',s['residual_evaluations']==budget)
                check(f'derivative_count_{name}_{budget}_{mode}',s['curvature_evaluations']==(budget-1)//2)
                check(f'capped_remainder_{name}_{budget}_{mode}',frac(s['remainder_upper'])<=frac(s['constant_remainder_upper']))
                consistent(name,s,f'fixed_reference_{name}_{budget}_{mode}')
                if previous is not None:
                    check(f'unchanged_z_{name}_{budget}_{mode}',s['z_lower']==previous['z_lower'] and s['z_upper']==previous['z_upper'])
                previous=s;screen.append(s)
        save('fixed_grid.json',screen)
        print('FIXED',name,flush=True)
    # Nonuniform nondyadic coordinates stress interval width handling.
    nonuniform=[]
    nodes=list(map(Fraction,['0','1/1000','1/100','1/2','3/5','1']))
    for name in ('sine_20','bump_0.005','PINN_k1_seed7'):
        for mode in ('constant','full','moment'):
            v=ShapedVerifier(fs[name](),mode);s=serialize(v.certify(nodes));s.update(problem=name)
            consistent(name,s,f'nonuniform_reference_{name}_{mode}');nonuniform.append(s)
    save('nonuniform.json',nonuniform)
    # Fixed final mesh size, all intermediate indicator calls counted.
    adaptive_screen=[]
    for name in ('sine_50','bump_0.005','PINN_k1_seed7'):
        for n in (128,256,512):
            for mode,indicator in [('constant','global'),('moment','global'),('moment','peak')]:
                v=ShapedVerifier(fs[name](),mode,indicator);v.refine_to(n);s=serialize(v.certify())
                s.update(problem=name,final_budget=n,mesh=[str(x) for x in v.nodes])
                consistent(name,s,f'adaptive_reference_{name}_{n}_{mode}_{indicator}')
                adaptive_screen.append(s)
        save('adaptive_grid.json',adaptive_screen)
    tasks=[('sine_50',Fraction(6,5)),('bump_0.005',Fraction(6,5)),
           ('bump_0.005',Fraction(51,50)),('PINN_k1_seed7',Fraction(501,500)),
           ('PINN_k5_seed7',Fraction(51,50))]
    methods=[('constant',False,'global'),('full',False,'global'),('moment',False,'global'),
             ('constant',True,'global'),('moment',True,'global'),('moment',True,'peak')]
    runs=[]
    for repeat in range(3):
        for name,target in tasks:
            ordered=methods[repeat:]+methods[:repeat]
            for mode,adaptive,indicator in ordered:
                s=request(fs[name],mode,adaptive,indicator,target);s.update(problem=name,repeat=repeat)
                consistent(name,s,f'timed_reference_{repeat}_{name}_{target}_{mode}_{adaptive}_{indicator}')
                if s['success']:check(f'target_{len(runs)}',frac(s['upper'])<=target*frac(s['lower']))
                runs.append(s);save('runs.json',runs)
                print('TIMED',repeat,name,str(target),mode,'adaptive',adaptive,indicator,
                      'N',s['intervals'],'ms',round(s['wall_ms'],2),'success',s['success'],flush=True)
    summary=[]
    for name,target in tasks:
        for mode,adaptive,indicator in methods:
            group=[r for r in runs if r['problem']==name and r['target']==str(target)
                   and (r['mode'],r['adaptive'],r['indicator'])==(mode,adaptive,indicator)]
            summary.append(dict(problem=name,target=str(target),mode=mode,adaptive=adaptive,indicator=indicator,
                                successes=sum(r['success'] for r in group),intervals=group[0]['intervals'],
                                residual_evaluations=group[0]['residual_evaluations'],curvature_evaluations=group[0]['curvature_evaluations'],
                                median_wall_ms=statistics.median(r['wall_ms'] for r in group),
                                min_wall_ms=min(r['wall_ms'] for r in group),max_wall_ms=max(r['wall_ms'] for r in group)))
    save('summary.json',summary)
    config.update(status='complete',checks_passed=len(checks),fixed_grid_configurations=len(screen),
                  nonuniform_configurations=len(nonuniform),adaptive_configurations=len(adaptive_screen),
                  timed_requests=len(runs),successful_requests=sum(r['success'] for r in runs))
    save('config.json',config)
    print('COMPLETED',out,'checks',len(checks),flush=True)


if __name__=='__main__':main()
