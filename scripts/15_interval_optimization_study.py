"""Controlled optimization study: interval expressions, native backend and mesh.

Uses frozen weights only. All stopping rules use certificate bounds, never the
known exact solution. The separate reference enclosures are diagnostic oracles.
No GitHub writes. Each run saves a new evidence directory and checkpoints.
"""
import argparse
from datetime import datetime,timezone
from fractions import Fraction
import gc
import json
import math
from pathlib import Path
import platform
import statistics
import sys
import time
import hashlib

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpmath as mp
import flint
from flint import arb,ctx as arb_ctx

from certification_core import NeuralResidual,certificate as float_certificate
from interval_certification import (context,FrozenNeuralPoisson,certify as mp_certify,
    endpoint_fraction,export_bounds,continuous_reference,sha256)
from interval_cached import CachedFrozenNeuralPoisson
from interval_reformulated import ReformulatedPoisson,factored_polynomials
from arb_certification import (ArbPoisson,certify as arb_certify,exact_fraction,
    bounds as arb_bounds,square,polys,export as arb_export)

ROOT=Path(__file__).resolve().parents[1]
MAIN=ROOT/'cases/poisson_1d_sine_k1_seed7'
AUDIT=ROOT/'results/2026-10-28/research_runs/12_interval_audit/20261028T125218065148Z'
VARIANTS=['mp_cached','mp_centered','mp_factored','mp_critical',
          'arb_expanded','arb_factored','arb_critical']
TARGETS=[Fraction(6,5),Fraction(51,50),Fraction(501,500)]


def frac(e):return Fraction(int(e['numerator']),int(e['denominator']))


def make(case,variant,bits=103):
    if variant.startswith('arb_'):return ArbPoisson(case,bits,variant.removeprefix('arb_'))
    c=context(30 if bits==103 else 50)
    if variant=='mp_cached':return CachedFrozenNeuralPoisson(case,c)
    return ReformulatedPoisson(case,c,variant.removeprefix('mp_'))


def solve(p,variant,n):
    return arb_certify(p,n) if variant.startswith('arb_') else mp_certify(p,n)


def serialize(r,variant):
    native=variant.startswith('arb_')
    result=arb_bounds(r) if native else export_bounds(r)
    point_export=arb_export if native else None
    if native:
        result['remainder_upper']=point_export(r['remainder_upper'],'up')
    else:
        from interval_certification import export_endpoint
        result['remainder_upper']=export_endpoint(r['remainder_upper'],'up')
    result['intervals']=r['intervals'];result['elapsed_s']=r['elapsed_s']
    result['gap_display']=float(frac(result['upper'])/frac(result['lower'])) if frac(result['lower'])>0 else None
    return result


def create_secondary_cases(out):
    fixture=json.loads((ROOT/'configs/frozen_networks_round07.json').read_text())
    template=json.loads((MAIN/'case.json').read_text())
    paths={'k1_seed7':MAIN}
    for seed in (17,27):
        row=next(r for r in fixture['networks'] if r['frequency']==1 and r['seed']==seed)
        directory=out/'fixtures'/f'k1_seed{seed}';directory.mkdir(parents=True)
        (directory/'weights.json').write_text(json.dumps(row['weights'],indent=2)+'\n')
        (directory/'provenance.json').write_text(json.dumps(dict(source='configs/frozen_networks_round07.json',
            source_sha256=sha256(ROOT/'configs/frozen_networks_round07.json'),frequency=1,seed=seed,
            training_steps=2000,reference_kind='historical observed grid maximum',
            reference_error_grid=row['reference_error_grid']),indent=2)+'\n')
        (directory/'identity.json').write_text(json.dumps(dict(label=row['label'],
            scope='Identity-only secondary case; no reused main-case certificate baseline'),indent=2)+'\n')
        c=json.loads(json.dumps(template));c['case_id']=f'poisson_1d_sine_k1_seed{seed}'
        c['network']['weights_sha256']=sha256(directory/'weights.json')
        c['provenance_file']='provenance.json';c['provenance_sha256']=sha256(directory/'provenance.json')
        c['baseline_file']='identity.json';c['baseline_sha256']=sha256(directory/'identity.json')
        c.pop('reference',None)
        c['selection_reason']='Additional frozen k=1 seed, no retraining; secondary consistency case'
        (directory/'case.json').write_text(json.dumps(c,indent=2)+'\n')
        paths[f'k1_seed{seed}']=directory
    return paths


def float_problem(case):
    c=json.loads((Path(case)/'case.json').read_text())
    w=json.loads((Path(case)/c['network']['weights_file']).read_text())
    return NeuralResidual((np.array(w['hidden.weight']).flatten(),np.array(w['hidden.bias']),
        np.array(w['output.weight']).flatten(),w['output.bias'][0]),1,c['case_id'])


def pilot(case,target):
    """Cheap heuristic only. No monotonicity assumption is needed for validity.

    Float gap bisection proposes a mesh with some margin. The interval result
    must independently pass; subsequent refinements always increase the mesh.
    """
    p=float_problem(case);goal=1+.75*(float(target)-1);lo=16;hi=32;queries=0
    def passes(n):
        nonlocal queries
        queries+=n
        r=float_certificate(p,np.linspace(0,1,n+1),'linear')
        return r['lower']>0 and r['upper']<=goal*r['lower']
    while not passes(hi):
        lo=hi;hi*=2
        if hi>4096:raise RuntimeError('Pilot cap exceeded')
    while hi-lo>8:
        mid=(hi+lo)//2
        if passes(mid):hi=mid
        else:lo=mid
    return max(32,8*math.ceil(hi/8)),queries


def timed_request(case,variant,controller,target):
    # Clear debt from the previous request outside the clock, but leave normal
    # collection enabled within this request and record its collection counts.
    gc.collect();before_gc=[s['collections'] for s in gc.get_stats()]
    start=time.perf_counter();cpu=time.process_time();float_cells=0
    n=32
    if controller=='pilot':n,float_cells=pilot(case,target)
    p=make(case,variant);history=[];count=0
    while n<=4096:
        r=solve(p,variant,n);s=serialize(r,variant);count+=n
        low,up=frac(s['lower']),frac(s['upper'])
        success=low>0 and up<=target*low
        history.append(s)
        if success:break
        if controller=='doubling':n*=2
        else:
            # A cost heuristic, not a certificate. Always verify the new mesh.
            gap=float(up/low) if low>0 else None
            factor=max(1.2,min(2.,1.08*math.sqrt((gap-1)/(float(target)-1)))) if gap else 2.
            n=8*math.ceil(n*factor/8)
    if not success:raise RuntimeError('Interval budget exhausted')
    result=dict(variant=variant,controller=controller,target=float(target),intervals=r['intervals'],
        lower=s['lower'],upper=s['upper'],gap_display=s['gap_display'],history=history,
        interval_cells=count,float_cells=float_cells,wall_ms=(time.perf_counter()-start)*1000,
        cpu_ms=(time.process_time()-cpu)*1000)
    result['gc_collections']=[s['collections']-b for s,b in zip(gc.get_stats(),before_gc)]
    return result


class AnalyticArb:
    boundary_error=(0,0)
    def __init__(self,name,bits=103):self.name=name;self.bits=bits;arb_ctx.prec=bits
    def guard(self):
        if arb_ctx.prec!=self.bits:raise ValueError('Precision changed')
    def residual(self,x):
        if self.name=='zero':return arb(0)
        if self.name=='constant':return arb(2)
        if self.name=='quartic':return (-32+192*x-192*square(x))/10000
        k=int(self.name.split('_')[1]);w=k*arb.pi()
        return w*w*(w*x).sin()/10000
    def r_second(self,x):
        if self.name in ('zero','constant'):return arb(0)
        if self.name=='quartic':return -arb(384)/10000
        k=int(self.name.split('_')[1]);w=k*arb.pi();w2=w*w
        return -w2*w2*(w*x).sin()/10000
    def target(self):
        if self.name=='zero':return Fraction(0)
        if self.name=='constant':return Fraction(1,4)
        return Fraction(1,10000)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path)
    parser.add_argument('--repeats',type=int,default=5);args=parser.parse_args()
    if args.repeats<3:raise ValueError('At least three repeats required')
    out=args.out or ROOT/'results/generated/15_interval_optimization'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True,exist_ok=False);checks=[]
    def check(name,ok,**details):
        checks.append(dict(name=name,passed=bool(ok),**details))
        (out/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
        if not ok:raise AssertionError(name)
    cases=create_secondary_cases(out);references={}
    references['k1_seed7']=json.loads((AUDIT/'continuous_reference.json').read_text())
    for key,case in cases.items():
        if key=='k1_seed7':continue
        r=continuous_reference(FrozenNeuralPoisson(case,context(50)),2048)
        references[key]=export_bounds(r)
        print('Continuous reference computed:',key,flush=True)
    (out/'continuous_references.json').write_text(json.dumps(references,indent=2)+'\n')
    # Exact symbolic coefficient identities, independent of floating arithmetic.
    def mul(a,b):
        c=[0]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):c[i+j]+=x*y
        return c
    exact_factors={2:mul([0,-2],[1,0,-1]),3:mul([-2,0,2],[1,0,-3]),
                   4:mul(mul([0,8],[1,0,-1]),[2,0,-3])}
    expected={2:[0,-2,0,2],3:[-2,0,8,0,-6],4:[0,16,0,-40,0,24]}
    for j in expected:check('exact_polynomial_identity_'+str(j),exact_factors[j]==expected[j])
    # Native conversion and crossing-zero square/serialization guards.
    for a,b in [(Fraction(-1),Fraction(1)),(Fraction(1,4),Fraction(3,4)),(Fraction(-3),Fraction(-2))]:
        arb_ctx.prec=103;x=(arb(a.numerator)/a.denominator).union(arb(b.numerator)/b.denominator);sq=square(x)
        true_lo=0 if a<=0<=b else min(a*a,b*b);true_hi=max(a*a,b*b)
        check('native_square_'+str(a)+'_'+str(b),exact_fraction(sq.lower())<=true_lo and exact_fraction(sq.upper())>=true_hi)
    check('native_unhandled_power_detected',not (arb(-1).union(1)**2).is_finite())
    class NonfiniteProblem:
        boundary_error=(0,0);bits=103
        def __init__(self,where,value):self.where=where;self.bad=arb(value)
        def guard(self):arb_ctx.prec=103
        def residual(self,x):return self.bad if self.where=='residual' else arb(0)
        def r_second(self,x):return self.bad if self.where=='curvature' else arb(0)
    for where in ('residual','curvature'):
        for bad in ('nan','+inf','-inf'):
            try:arb_certify(NonfiniteProblem(where,bad),4)
            except (ValueError,ArithmeticError):rejected=True
            else:rejected=False
            check('reject_nonfinite_'+where+'_'+bad,rejected)
    guarded=ArbPoisson(MAIN,103,'factored');arb_ctx.prec=169
    try:guarded.residual(arb(1)/2)
    except ValueError:rejected=True
    else:rejected=False
    check('native_precision_change_rejected',rejected)
    arb_ctx.prec=103
    for value in [arb(0),(arb(1)/3).lower(),-(arb(1)/7).upper(),arb('1e-100').lower()]:
        for direction in ('down','up'):
            e=arb_export(value,direction);q=exact_fraction(value);d=Fraction(e['decimal_outward'])
            check('native_export_'+str(len(checks)),frac(e)==q and (d<=q if direction=='down' else d>=q))
    for name,ns in [('zero',[1,17]),('constant',[1,17,64]),('quartic',[17,128]),('sine_1',[64]),('sine_20',[256]),('sine_200',[100])]:
        p=AnalyticArb(name)
        for n in ns:
            r=arb_certify(p,n);check(f'known_norm_{name}_{n}',exact_fraction(r['lower'])<=p.target()<=exact_fraction(r['upper']))
    print('Symbolic identities, native square and known-norm checks passed.',flush=True)

    # 100-digit point and whole-cell derivative diagnostics for all three seeds.
    with mp.workdps(100):
        for key,case in cases.items():
            raw=json.loads((case/'weights.json').read_text())
            w=[mp.mpf(float(v[0])) for v in raw['hidden.weight']];b=[mp.mpf(float(v)) for v in raw['hidden.bias']]
            a=[mp.mpf(float(v)) for v in raw['output.weight'][0]];c=mp.mpf(float(raw['output.bias'][0]))
            def value(x):return x*(1-x)*(c+sum(aa*mp.tanh(ww*x+bb) for ww,bb,aa in zip(w,b,a)))
            for variant in ['mp_factored','mp_critical','arb_factored','arb_critical']:
                p=make(case,variant)
                for numerator,denominator in [(0,1),(1,10),(1,2),(9,10),(1,1)]:
                    xx=mp.mpf(numerator)/denominator
                    xi=arb(numerator)/denominator if variant.startswith('arb_') else p.ctx.mpf(numerator)/denominator
                    vals=[('value',p.solution(xi),value(xx)),('r',p.residual(xi),-mp.diff(value,xx,2)-mp.pi**2*mp.sin(mp.pi*xx)),
                          ('r2',p.r_second(xi),-mp.diff(value,xx,4)+mp.pi**4*mp.sin(mp.pi*xx))]
                    for label,enclosure,v in vals:
                        q=Fraction(*mp.libmp.to_rational(v._mpf_))
                        lo=exact_fraction(enclosure.lower()) if variant.startswith('arb_') else endpoint_fraction(enclosure)
                        hi=exact_fraction(enclosure.upper()) if variant.startswith('arb_') else endpoint_fraction(enclosure,1)
                        check(f'point_{key}_{variant}_{numerator}_{denominator}_{label}',lo<=q<=hi)
                for i in (0,3,7):
                    if variant.startswith('arb_'):cell=(arb(i)/8).union(arb(i+1)/8)
                    else:cell=p.ctx.mpf([p.ctx.mpf(i)/8,p.ctx.mpf(i+1)/8])
                    enclosure=p.r_second(cell)
                    for j in (0,1,2):
                        xx=mp.mpf(2*i+j)/16;v=-mp.diff(value,xx,4)+mp.pi**4*mp.sin(mp.pi*xx)
                        q=Fraction(*mp.libmp.to_rational(v._mpf_))
                        lo=exact_fraction(enclosure.lower()) if variant.startswith('arb_') else endpoint_fraction(enclosure)
                        hi=exact_fraction(enclosure.upper()) if variant.startswith('arb_') else endpoint_fraction(enclosure,1)
                        check(f'cell_{key}_{variant}_{i}_{j}',lo<=q<=hi)
    print('High-precision diagnostics passed for all three frozen seeds.',flush=True)

    screening=[]
    for variant in VARIANTS:
        p=make(MAIN,variant)
        for n in (64,128,256,512):
            r=solve(p,variant,n);s=serialize(r,variant);s['variant']=variant
            ref=references['k1_seed7'];check(f'screen_reference_{variant}_{n}',frac(s['lower'])<=frac(ref['lower']) and frac(s['upper'])>=frac(ref['upper']))
            screening.append(s)
        (out/'screening.json').write_text(json.dumps(screening,indent=2)+'\n')
        print('Screened',variant,'512-cell gap',screening[-1]['gap_display'],'seconds',screening[-1]['elapsed_s'],flush=True)
    secondary=[]
    for key,case in cases.items():
        for variant in ['mp_cached','mp_factored','arb_factored','arb_critical']:
            if key=='k1_seed7':continue
            r=solve(make(case,variant),variant,256);s=serialize(r,variant);s.update(case=key,variant=variant)
            ref=references[key];check(f'secondary_reference_{key}_{variant}',frac(s['lower'])<=frac(ref['lower']) and frac(s['upper'])>=frac(ref['upper']))
            secondary.append(s)
    (out/'secondary_cases.json').write_text(json.dumps(secondary,indent=2)+'\n')
    # Additional precision check for the new native critical range backend.
    for key,case in cases.items():
        r=solve(make(case,'arb_critical',169),'arb_critical',256);s=serialize(r,'arb_critical')
        ref=references[key];check('native_169bit_reference_'+key,frac(s['lower'])<=frac(ref['lower']) and frac(s['upper'])>=frac(ref['upper']))
    print('Fixed-budget and secondary-case checks passed. Starting paired timing.',flush=True)

    methods=[('mp_cached','doubling'),('mp_factored','doubling'),('arb_factored','doubling'),
             ('arb_critical','doubling'),('mp_cached','pilot'),('arb_factored','pilot'),('arb_critical','pilot')]
    # Warm code paths, with no large context debt carried into timed requests.
    for variant,controller in methods:timed_request(MAIN,variant,controller,TARGETS[0])
    runs=[]
    for repeat in range(args.repeats):
        for target in TARGETS:
            rotated=methods[repeat%len(methods):]+methods[:repeat%len(methods)]
            for variant,controller in rotated:
                r=timed_request(MAIN,variant,controller,target);r['repeat']=repeat
                ref=references['k1_seed7']
                check(f'timed_ref_{repeat}_{target}_{variant}_{controller}',frac(r['lower'])<=frac(ref['lower']) and frac(r['upper'])>=frac(ref['upper']))
                check(f'timed_stop_{repeat}_{target}_{variant}_{controller}',frac(r['lower'])>0 and frac(r['upper'])<=target*frac(r['lower']))
                runs.append(r);(out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n')
                print(f"repeat {repeat+1} target {float(target)} {variant}/{controller}: n={r['intervals']} wall={r['wall_ms']:.2f}ms cpu={r['cpu_ms']:.2f}ms",flush=True)
    summary=[]
    for target in TARGETS:
        for variant,controller in methods:
            rows=[r for r in runs if r['target']==float(target) and r['variant']==variant and r['controller']==controller]
            summary.append(dict(target=float(target),variant=variant,controller=controller,intervals=rows[0]['intervals'],
                interval_cells=rows[0]['interval_cells'],float_cells=rows[0]['float_cells'],gap_display=rows[0]['gap_display'],
                median_wall_ms=statistics.median(r['wall_ms'] for r in rows),min_wall_ms=min(r['wall_ms'] for r in rows),
                max_wall_ms=max(r['wall_ms'] for r in rows),median_cpu_ms=statistics.median(r['cpu_ms'] for r in rows)))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(12,4.8))
    for variant in VARIANTS:
        rows=[r for r in screening if r['variant']==variant]
        axes[0].loglog([r['intervals'] for r in rows],[r['gap_display']-1 if r['gap_display'] else float('nan') for r in rows],'-o',label=variant,ms=3)
    axes[0].set(xlabel='Uniform intervals',ylabel='Upper / lower - 1',title='Same fixed weights: enclosure tightness');axes[0].legend(fontsize=6)
    labels=['MP cached','MP factored','Arb factored','Arb critical','MP pilot','Arb factored pilot','Arb critical pilot']
    x=np.arange(len(TARGETS));width=.11
    for i,((variant,controller),label) in enumerate(zip(methods,labels)):
        rows=[next(r for r in summary if r['target']==float(t) and r['variant']==variant and r['controller']==controller) for t in TARGETS]
        y=np.array([r['median_wall_ms'] for r in rows]);err=np.array([[r['median_wall_ms']-r['min_wall_ms'] for r in rows],[r['max_wall_ms']-r['median_wall_ms'] for r in rows]])
        axes[1].bar(x+(i-3)*width,y,width,yerr=err,capsize=2,label=label)
    axes[1].set(yscale='log',xticks=x,xticklabels=[str(float(t)) for t in TARGETS],xlabel='Required upper / lower',ylabel='Total wall time (ms)',title='Median and min/max; setup and attempts included');axes[1].legend(fontsize=6)
    for ax in axes:ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(out/'optimization.png',dpi=160);plt.close(fig)
    sources=[ROOT/'scripts'/n for n in ['15_interval_optimization_study.py','interval_reformulated.py','arb_certification.py','interval_certification.py','interval_cached.py','certification_core.py']]
    metadata=dict(date='2026-10-29',repeats=args.repeats,mpmath=mp.__version__,python_flint=flint.__version__,flint=flint.__FLINT_VERSION__,
        python=sys.version,platform=platform.platform(),gc_thresholds=gc.get_threshold(),checks_passed=len(checks),screening_configurations=len(screening),
        timed_requests=len(runs),sources={str(p.relative_to(ROOT)):sha256(p) for p in sources},case_sha256=sha256(MAIN/'case.json'),
        timing='Sequential rotated methods, fresh input and setup inside timer; gc.collect before each request outside timer; normal GC remains enabled and counts recorded. No exact solution used for mesh choice or stopping.',
        precision='mpmath 30 digits (103 endpoint bits), Arb 103 midpoint bits and its documented 30-bit radius; 169-bit native spot checks. Same acceptance criterion, not identical arithmetic formats.',
        limitations=['Three frozen k=1 shallow networks; no claim of all PINNs/PDEs','Assistant tests, no independent/formal software audit',
                     'Backend correctness and reviewed math remain assumptions','Native backend deliberately avoids zero-crossing arb ** integer observed indeterminate; uses enclosing square/multiplication',
                     'Do not compare yesterday timing values directly; GC/run conditions and machine load differ'])
    (out/'config.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print('Completed checks:',len(checks),'timed requests:',len(runs),flush=True)
    print('RESULT_DIRECTORY',out,flush=True)


if __name__=='__main__':main()
