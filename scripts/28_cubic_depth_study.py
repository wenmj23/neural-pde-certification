"""Cubic residual follow-up: equal samples, independent diagnostics, full costs."""
from datetime import datetime,timezone
from fractions import Fraction as F
import gc
import hashlib
import importlib.util
import json
from pathlib import Path
import statistics
import time
import mpmath as mp
from flint import arb

from arb_cubic_reconstruction import (CubicVerifier,SixthOrderDeep,SineSix,BumpSix,
    LAGRANGE,NODES,PRODUCT,CUBIC_AVERAGE,CUBIC_CORRECTION)
from arb_shaped_remainder import ShapedVerifier
from arb_problem_families import SmoothProblem,QuarticError
from arb_certification import exact_fraction,export

ROOT=Path(__file__).resolve().parents[1]
FIXTURE=ROOT/'configs/frozen_networks_depth_2026-10-29.json'
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
base=load('study25',ROOT/'scripts/25_depth_certification_study.py')
checks23=load('checks23',ROOT/'scripts/23_depth_derivative_checks.py')


def evaluate(poly,x):
    result=poly[0]*0
    for v in reversed(poly):result=result*x+v
    return result


class QuarticSix(QuarticError):
    def r_fourth(self,x):self.guard();return arb(0)


class CubicSource(SmoothProblem):
    def __init__(self,zero=False):self.setup(103);self.zero=zero
    def residual(self,x):self.guard();return arb(0) if self.zero else x*x*x
    def r_fourth(self,x):self.guard();return arb(0)


def request(factory,method,target,cap=1024,seconds=8):
    gc.collect();start=time.perf_counter();cpu=time.process_time();p=factory()
    v=ShapedVerifier(p,'moment') if method=='quadratic_moment' else CubicVerifier(p,'constant' if method=='cubic_constant' else 'moment')
    degree=2 if method=='quadratic_moment' else 3;m=4;history=[];last=None;success=False;status='mesh_cap'
    while m<=cap:
        v.uniform(6*m//degree)
        try:
            last=base.serialize(v.certify());lo,hi=base.frac(last['lower']),base.frac(last['upper']);success=lo>0 and hi<=target*lo
            history.append(dict(sample_budget=6*m+1,status='finite',lower=last['lower'],upper=last['upper']))
        except ArithmeticError as error:history.append(dict(sample_budget=6*m+1,status='safe_rejection',reason=str(error)))
        if success:status='success';break
        if time.perf_counter()-start>=seconds:status='time_cap';break
        m*=2
    return dict(method=method,target=str(target),success=bool(success),status=status,final=last,history=history,
                wall_ms=(time.perf_counter()-start)*1000,cpu_ms=(time.process_time()-cpu)*1000)


def main():
    out=ROOT/'results/generated/28_cubic_depth'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');out.mkdir(parents=True,exist_ok=False)
    sha=hashlib.sha256(FIXTURE.read_bytes()).hexdigest();rows=json.loads(FIXTURE.read_text())['networks']
    ref_run=max(p.parent for p in (ROOT/'results/generated/24_depth_reference').glob('*/config.json') if json.loads(p.read_text()).get('status')=='complete')
    refs=json.loads((ref_run/'references.json').read_text());checks=[];mp.mp.dps=180
    sources=[Path(__file__),ROOT/'scripts/arb_cubic_reconstruction.py',ROOT/'scripts/arb_deep_network.py',ROOT/'scripts/arb_shaped_remainder.py',ROOT/'scripts/arb_reconstruction.py',ROOT/'scripts/25_depth_certification_study.py',ROOT/'scripts/23_depth_derivative_checks.py',FIXTURE,ref_run/'references.json']
    (out/'source_snapshot').mkdir()
    for p in sources:
        if p.suffix=='.py':(out/'source_snapshot'/p.name).write_bytes(p.read_bytes())
    config=dict(status='running',bits=103,reference_decimal_digits=180,repeats=3,shape_average=str(CUBIC_AVERAGE),shape_correction=str(CUBIC_CORRECTION),
                protocol='Equal sample B=6m+1: quadratic 3m cells, cubic 2m cells; identical sample-budget doubling search m=4,8,...,1024, 8s soft cap checked between meshes.',
                regularity='Cubic requires residual C4, e.g. candidate C6 and forcing C4; tanh and sine fixtures satisfy this.',
                timing='Fresh model/cache and all budget attempts included; known error/reference checks and file writes excluded.',
                scope='Experimental higher order; conditional interval correctness, not independently/formally audited',
                sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    def save(name,value):(out/name).write_text(json.dumps(value,indent=2)+'\n')
    save('config.json',config)
    def check(name,passed):
        checks.append(dict(name=name,passed=bool(passed)));save('checks.json',checks)
        if not passed:raise AssertionError(name)
    check('product_center_exact',evaluate(PRODUCT,F(1,2))==F(1,144))
    q=lambda w:w*w-F(5,18)*w+F(1,144)
    check('product_absolute_maximum_exact',q(F(5,36))==-F(1,81) and q(F(0))==F(1,144) and q(F(1,4))==0)
    check('average_exact',CUBIC_AVERAGE==F(49,90))
    for i,poly in enumerate(LAGRANGE):
        for j,node in enumerate(NODES):check(f'lagrange_{i}_{j}',evaluate(poly,node)==int(i==j))
    # Exact cubic residual: interpolation remainder is zero, true error norm known.
    for mode in ('constant','moment'):
        for n in (1,4,16):
            v=CubicVerifier(CubicSource(),mode);v.uniform(n);r=v.certify();ref=1/(25*arb(5).sqrt().sqrt())
            check(f'exact_cubic_reconstruction_{mode}_{n}',r['lower']<=ref.upper() and r['upper']>=ref.lower() and exact_fraction(r['remainder_upper'])==0)
        v=CubicVerifier(CubicSource(True),mode);r=v.certify();check('zero_'+mode,exact_fraction(r['lower'])==exact_fraction(r['upper'])==0)
        v=CubicVerifier(CubicSource(),mode);r=v.certify([F(0),F(1,7),F(1,2),F(3,4),F(1)])
        ref=1/(25*arb(5).sqrt().sqrt());check('nonuniform_exact_cubic_'+mode,r['lower']<=ref.upper() and r['upper']>=ref.lower())
    alias=CubicVerifier(SineSix(120));alias.uniform(40);r=alias.certify()
    check('alias_samples_zero',max(exact_fraction(v.abs_upper()) for v in alias.node_cache.values())<F(1,10**20))
    check('alias_not_missed',exact_fraction(r['upper'])>=F(1,10000))
    class Invalid(CubicSource):
        def r_fourth(self,x):return arb('nan')
    rejected=False
    try:CubicVerifier(Invalid()).certify()
    except ArithmeticError:rejected=True
    check('nonfinite_rejected',rejected)
    # Validate the newly used sixth network derivative and r'''' at points.
    for row in rows:
        function=checks23.mp_function(row);omega=row['frequency']*mp.pi
        p=SixthOrderDeep(FIXTURE,row['label'],sha)
        for point in (F(1,7),F(1,2),F(6,7)):
            x=mp.mpf(point.numerator)/point.denominator;co=mp.taylor(function,x,6);xx=arb(point.numerator)/point.denominator;jet=p.network_jet(xx,6)
            for j,value in enumerate(jet):check(f'sixth_network_{row["label"]}_{point}_{j}',value.overlaps(arb(str(co[j]))))
            # v^(6)=720*(q*n6+q'*n5-n4) for normalized Taylor coefficients.
            reference=-720*(x*(1-x)*co[6]+(1-2*x)*co[5]-co[4])-omega**6*mp.sin(omega*x)
            check(f'fourth_residual_{row["label"]}_{point}',p.r_fourth(xx).overlaps(arb(str(reference))))
            if point==F(1,7):
                with mp.workdps(240):
                    hi_function=checks23.mp_function(row);hi=mp.taylor(hi_function,mp.mpf(1)/7,6)[6]
                    scale=max(abs(hi),abs(co[6]))
                    check(f'sixth_reference_convergence_{row["label"]}',abs(hi-co[6])<=scale*mp.mpf('1e-80') if scale else hi==co[6])
        print('SIXTH_CHECKED',row['label'],flush=True)
    factories={f'sine_{k}':(lambda k=k:SineSix(k)) for k in (1,20,50)}
    factories['mixed']=lambda:SineSix(1,frequencies=[1,5,81],coefficients=[F(120,128),F(7,128),F(1,128)])
    factories.update(bump_002=lambda:BumpSix('0.02'),bump_0005=lambda:BumpSix('0.005'),quartic=QuarticSix)
    for label,factory in factories.items():
        p=factory();function=p.mp_function()
        points=[F(1,7),F(2,5)] if not label.startswith('bump') else [F(1,2)+p.width/5,F(1,2)+2*p.width/5]
        for point in points:
            x=mp.mpf(point.numerator)/point.denominator
            reference=mp.mpf(0) if label=='quartic' else -mp.diff(function,x,6)
            check(f'analytic_sixth_{label}_{point}',p.r_fourth(arb(point.numerator)/point.denominator).overlaps(arb(str(reference))))
    for row in rows:
        label=row['label'];factories[label]=lambda label=label:SixthOrderDeep(FIXTURE,label,sha)
    def consistent(label,s,name):
        if label in refs:lo,hi=base.frac(refs[label]['lower']),base.frac(refs[label]['upper'])
        else:lo=hi=F(1,10000)
        check(name,base.frac(s['lower'])<=lo and base.frac(s['upper'])>=hi)
    methods=['quadratic_moment','cubic_constant','cubic_moment'];fixed=[]
    for label,factory in factories.items():
        for m in (8,32,128):
            for method in methods:
                p=factory();degree=2 if method=='quadratic_moment' else 3
                v=ShapedVerifier(p,'moment') if degree==2 else CubicVerifier(p,'constant' if method=='cubic_constant' else 'moment')
                v.uniform(6*m//degree);s=base.serialize(v.certify());consistent(label,s,f'fixed_ref_{label}_{m}_{method}')
                check(f'equal_samples_{label}_{m}_{method}',s['residual_evaluations']==6*m+1)
                fixed.append(dict(label=label,method=method,sample_budget=6*m+1,**s))
        save('fixed_grid.json',fixed);print('CUBIC_FIXED',label,flush=True)
    tasks=[(r['label'],F(51,50)) for r in rows]+[('sine_50',F(6,5)),('bump_0005',F(6,5)),('bump_0005',F(51,50))]
    runs=[]
    for repeat in range(3):
        for label,target in tasks:
            for method in methods[repeat:]+methods[:repeat]:
                r=request(factories[label],method,target);r.update(label=label,repeat=repeat)
                if r['final'] is not None:consistent(label,r['final'],f'timed_ref_{repeat}_{label}_{target}_{method}')
                runs.append(r);save('runs.json',runs)
                print('CUBIC_TIMED',repeat,label,str(target),method,r['status'],'N',r['final']['intervals'] if r['final'] else None,'ms',round(r['wall_ms'],1),flush=True)
    summary=[]
    for label,target in tasks:
        for method in methods:
            group=[r for r in runs if r['label']==label and r['target']==str(target) and r['method']==method]
            summary.append(dict(label=label,target=str(target),method=method,successes=sum(r['success'] for r in group),
                                intervals=group[0]['final']['intervals'] if group[0]['final'] else None,
                                median_wall_ms=statistics.median(r['wall_ms'] for r in group),min_wall_ms=min(r['wall_ms'] for r in group),max_wall_ms=max(r['wall_ms'] for r in group)))
    save('summary.json',summary);config.update(status='complete',checks_passed=len(checks),fixed_configurations=len(fixed),timed_requests=len(runs),successful_requests=sum(r['success'] for r in runs));save('config.json',config)
    print('COMPLETED',out,'checks',len(checks),'success',config['successful_requests'],'/',len(runs),flush=True)


if __name__=='__main__':main()
