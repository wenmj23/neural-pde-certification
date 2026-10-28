"""Matched-gap costs of float, interval, caching and float-guided verification.

All timed modes reload the same frozen case; no training or exact solution is
used to select mesh sizes or stop. Includes setup, failed refinements, and safe
interval endpoint export. Excludes imports, result-file I/O and research-only
reference checks. Sequential/interleaved repeats, wall and process CPU clocks.
"""
import argparse
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import statistics
import sys
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpmath
from scipy.linalg import solve_banded
import scipy

from certification_core import NeuralResidual, certificate
from interval_certification import (context, FrozenNeuralPoisson, certify,
    endpoint_fraction, export_bounds, signed_polynomial_norm,
    remainder_norm_upper, mesh, sha256)
from interval_cached import CachedFrozenNeuralPoisson

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/'cases/poisson_1d_sine_k1_seed7'
AUDIT = ROOT/'results/generated/12_interval_audit/20261028T125218065148Z'
TARGETS = [Fraction(6,5), Fraction(51,50), Fraction(501,500)]
MODES = ['float64', 'interval_original', 'interval_cached', 'float_pilot_interval_cached']


def fraction(e):
    return Fraction(int(e['numerator']),int(e['denominator']))


def load_float():
    case=json.loads((CASE/'case.json').read_text())
    for filename,digest in [(case['network']['weights_file'],case['network']['weights_sha256']),
                            (case['provenance_file'],case['provenance_sha256']),
                            (case['baseline_file'],case['baseline_sha256'])]:
        if sha256(CASE/filename)!=digest: raise ValueError('Frozen case changed')
    if case['problem']['frequency']!=1 or case['network']['architecture']!=[1,16,1]:
        raise ValueError('Unexpected fixed model')
    w=json.loads((CASE/case['network']['weights_file']).read_text())
    params=(np.array(w['hidden.weight']).flatten(),np.array(w['hidden.bias']),
            np.array(w['output.weight']).flatten(),float(w['output.bias'][0]))
    if [np.size(x) for x in params]!=[16,16,16,1] or not all(np.isfinite(x).all() for x in params):
        raise ValueError('Invalid fixed weights')
    return NeuralResidual(params,1,case['case_id'])


def float_search(problem,target,cap=2048):
    n=32;history=[]
    while n<=cap:
        r=certificate(problem,np.linspace(0,1,n+1),'linear')
        success=r['lower']>0 and r['upper']<=float(target)*r['lower']
        history.append(dict(stage='float',intervals=n,lower=r['lower'],upper=r['upper']))
        if success or n==cap: break
        n*=2
    return r,history,bool(success)


def interval_search(problem,target,start=32,cap=2048):
    n=start;history=[]
    while n<=cap:
        r=certify(problem,n)
        low,up=endpoint_fraction(r['lower']),endpoint_fraction(r['upper'])
        success=low>0 and up<=target*low
        history.append(dict(stage='interval',intervals=n,**export_bounds(r)))
        if success or n==cap: break
        n*=2
    return r,history,bool(success)


def run(mode,target):
    wall=time.perf_counter();cpu=time.process_time();history=[];setup_ms=0
    if mode=='float64':
        before=time.perf_counter();p=load_float();setup_ms+=(time.perf_counter()-before)*1000
        r,history,success=float_search(p,target)
        bounds=dict(lower=r['lower'],upper=r['upper'])
        n=history[-1]['intervals']
        gap=r['upper']/r['lower'] if r['lower']>0 else None
    else:
        n=32
        if mode=='float_pilot_interval_cached':
            before=time.perf_counter();p=load_float();setup_ms+=(time.perf_counter()-before)*1000
            _,pilot,pilot_success=float_search(p,target)
            history.extend(pilot);n=pilot[-1]['intervals']
            if not pilot_success: raise RuntimeError('Float pilot reached its cap')
        before=time.perf_counter()
        cls=FrozenNeuralPoisson if mode=='interval_original' else CachedFrozenNeuralPoisson
        p=cls(CASE,context(30));setup_ms+=(time.perf_counter()-before)*1000
        r,trials,success=interval_search(p,target,start=n)
        history.extend(trials)
        bounds=export_bounds(r);n=trials[-1]['intervals']
        gap=float(endpoint_fraction(r['upper'])/endpoint_fraction(r['lower'])) if endpoint_fraction(r['lower'])>0 else None
    wall_ms=(time.perf_counter()-wall)*1000
    cpu_ms=(time.process_time()-cpu)*1000
    return dict(mode=mode,target_numerator=target.numerator,target_denominator=target.denominator,
                target=float(target),success=success,intervals=n,gap_display=gap,
                **bounds,wall_ms=wall_ms,cpu_ms=cpu_ms,setup_ms=setup_ms,history=history,
                float_cells=sum(t['intervals'] for t in history if t['stage']=='float'),
                interval_cells=sum(t['intervals'] for t in history if t['stage']=='interval'),
                conditional_interval_certificate=mode!='float64')


def profile_interval(cls,intervals=512):
    start=time.perf_counter();cpu=time.process_time();p=cls(CASE,context(30))
    times={'setup':(time.perf_counter()-start)*1000};ctx=p.ctx
    before=time.perf_counter();nodes=mesh(ctx,intervals);times['mesh']=(time.perf_counter()-before)*1000
    before=time.perf_counter();r=[p.residual(x) for x in nodes];times['point_residual']=(time.perf_counter()-before)*1000
    before=time.perf_counter();delta=[]
    for a,b in zip(nodes[:-1],nodes[1:]):
        K=abs(p.r_second(ctx.mpf([a.a,b.b]))).b
        delta.append((K*(b-a)**2/8).b)
    times['cell_second_derivative']=(time.perf_counter()-before)*1000
    before=time.perf_counter();zl,zu=signed_polynomial_norm(ctx,nodes,r);times['signed_integral']=(time.perf_counter()-before)*1000
    before=time.perf_counter();du=remainder_norm_upper(ctx,nodes,delta);times['remainder_integral']=(time.perf_counter()-before)*1000
    before=time.perf_counter();result=dict(lower=max(ctx.mpf(0),(zl-du).a),upper=(zu+du).b)
    bounds=export_bounds(result);times['final_export']=(time.perf_counter()-before)*1000
    return dict(implementation=cls.__name__,intervals=intervals,times_ms=times,
                wall_ms=(time.perf_counter()-start)*1000,cpu_ms=(time.process_time()-cpu)*1000,**bounds)


def finite_difference(intervals):
    h=1/intervals;x=np.linspace(0,1,intervals+1)
    band=np.zeros((3,intervals-1));band[0,1:]=-1;band[1,:]=2;band[2,:-1]=-1
    values=solve_banded((1,1),band,h*h*np.pi**2*np.sin(np.pi*x[1:-1]))
    return x,np.r_[0.,values,0.]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--repeats',type=int,default=3);parser.add_argument('--out',type=Path)
    args=parser.parse_args()
    if args.repeats<3: raise ValueError('Use at least three repeats for this preliminary benchmark')
    out=args.out or ROOT/'results/generated/13_verification_cost' / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True,exist_ok=False)
    reference=json.loads((AUDIT/'continuous_reference.json').read_text())
    audited=json.loads((AUDIT/'certificates.json').read_text())
    checks=[]
    def check(name,condition):
        checks.append(dict(name=name,passed=bool(condition)))
        (out/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
        if not condition: raise AssertionError(name)
    for dps in (30,50):
        p=CachedFrozenNeuralPoisson(CASE,context(dps))
        for n in (128,512,1000):
            new=certify(p,n)
            previous=next(r for r in audited if r['precision_digits']==dps and r['intervals']==n)
            for key,oldkey in [('lower','lower'),('upper','upper'),('z_lower','auxiliary_norm_lower'),
                               ('z_upper','auxiliary_norm_upper'),('remainder_upper','remainder_norm_upper')]:
                check(f'cached_exact_identity_{dps}_{n}_{key}',endpoint_fraction(new[key])==fraction(previous[oldkey]))
            print('Cached endpoint identity checked:',dps,'digits,',n,'intervals',flush=True)
    p=CachedFrozenNeuralPoisson(CASE,context(30));p.ctx.dps=40
    try:p.residual(p.ctx.mpf(1)/2)
    except ValueError:guard=True
    else:guard=False
    check('precision_mutation_rejected',guard)
    fallback=run('float_pilot_interval_cached',Fraction(23,20))
    attempts=[r for r in fallback['history'] if r['stage']=='interval']
    check('pilot_is_proposal_final_interval_must_pass',fallback['success'] and len(attempts)>=2
          and attempts[0]['intervals']==128 and fallback['intervals']==256)
    (out/'pilot_fallback_check.json').write_text(json.dumps(fallback,indent=2)+'\n')
    print('Fallback test passed: interval check rejects 128-cell pilot, then refines to 256.',flush=True)
    # Warm only implementations/code paths; main runs always start with fresh
    # candidate objects. Correctness runs above are excluded from timings.
    for mode in MODES:run(mode,Fraction(6,5))
    print('Timed benchmark begins; all modes run sequentially.',flush=True)
    records=[]
    for repeat in range(args.repeats):
        for target in TARGETS:
            order=MODES[repeat%len(MODES):]+MODES[:repeat%len(MODES)]
            for mode in order:
                r=run(mode,target);r['repeat']=repeat
                check(f'run_success_{repeat}_{target}_{mode}',r['success'])
                if r['conditional_interval_certificate']:
                    check(f'continuous_reference_contained_{repeat}_{target}_{mode}',
                          fraction(r['lower'])<=fraction(reference['lower'])
                          and fraction(r['upper'])>=fraction(reference['upper']))
                    check(f'exact_gap_stop_{repeat}_{target}_{mode}',
                          fraction(r['lower'])>0 and fraction(r['upper'])<=target*fraction(r['lower']))
                records.append(r)
                (out/'runs.json').write_text(json.dumps(records,indent=2)+'\n')
                print(f"repeat={repeat+1} target={float(target)} {mode}: cells={r['intervals']} wall={r['wall_ms']:.3f}ms cpu={r['cpu_ms']:.3f}ms",flush=True)
    summary=[]
    for target in TARGETS:
        for mode in MODES:
            rs=[r for r in records if r['mode']==mode and r['target']==float(target)]
            check(f'deterministic_grid_{target}_{mode}',len({r['intervals'] for r in rs})==1)
            summary.append(dict(mode=mode,target=float(target),intervals=rs[0]['intervals'],
                                gap_display=rs[0]['gap_display'],float_cells=rs[0]['float_cells'],
                                interval_cells=rs[0]['interval_cells'],
                                median_wall_ms=statistics.median(r['wall_ms'] for r in rs),
                                min_wall_ms=min(r['wall_ms'] for r in rs),max_wall_ms=max(r['wall_ms'] for r in rs),
                                median_cpu_ms=statistics.median(r['cpu_ms'] for r in rs),
                                conditional_interval_certificate=mode!='float64'))
        iv=[r for r in records if r['target']==float(target) and r['conditional_interval_certificate']]
        check(f'all_interval_modes_same_exact_endpoints_{target}',len({(fraction(r['lower']),fraction(r['upper'])) for r in iv})==1)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    profiles=[]
    for _ in range(3):
        for cls in (FrozenNeuralPoisson,CachedFrozenNeuralPoisson):
            r=profile_interval(cls)
            prior=next(p for p in audited if p['precision_digits']==30 and p['intervals']==512)
            check('profile_endpoint_identity_'+str(len(profiles)),fraction(r['lower'])==fraction(prior['lower'])
                  and fraction(r['upper'])==fraction(prior['upper']))
            profiles.append(r)
    (out/'profiles.json').write_text(json.dumps(profiles,indent=2)+'\n')
    fd=[]
    for n in (100,1000):
        for _ in range(3):finite_difference(n)
        samples=[]
        for _ in range(30):
            w=time.perf_counter();c=time.process_time();finite_difference(n)
            samples.append(dict(wall_ms=(time.perf_counter()-w)*1000,cpu_ms=(time.process_time()-c)*1000))
        x,y=finite_difference(n);xx=np.linspace(0,1,10001)
        fd.append(dict(intervals=n,median_wall_ms=statistics.median(s['wall_ms'] for s in samples),
                       median_cpu_ms=statistics.median(s['cpu_ms'] for s in samples),samples=samples,
                       observed_interpolated_grid_error=float(np.max(abs(np.interp(xx,x,y)-np.sin(np.pi*xx)))),
                       caveat='Solves the PDE; does not certify the given PINN. Ordinary floating result.'))
    (out/'finite_difference.json').write_text(json.dumps(fd,indent=2)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(12,4.8))
    labels=['Float64 (diagnostic)','Interval original','Interval cached','Float pilot + cached interval']
    colors=['#577590','#f3722c','#43aa8b','#9b5de5'];positions=np.arange(len(TARGETS));width=.2
    for j,(mode,label,color) in enumerate(zip(MODES,labels,colors)):
        vals=[next(r['median_wall_ms'] for r in summary if r['mode']==mode and r['target']==float(t)) for t in TARGETS]
        axes[0].bar(positions+(j-1.5)*width,vals,width,label=label,color=color)
    axes[0].set(yscale='log',ylabel='Median total wall time (ms, log scale)',xlabel='Required upper / lower ratio',
                title='Same frozen PINN, refinement attempts included',xticks=positions,xticklabels=[str(float(t)) for t in TARGETS])
    axes[0].legend(fontsize=7)
    stages=['setup','mesh','point_residual','cell_second_derivative','signed_integral','remainder_integral','final_export']
    stage_labels=['Setup','Mesh','Point residuals','Cell residual curvature','Signed integral','Remainder integral','Export']
    stage_colors=['#adb5bd','#c9ada7','#577590','#f3722c','#43aa8b','#9b5de5','#f9c74f']
    bottom=np.zeros(2)
    for stage,label,color in zip(stages,stage_labels,stage_colors):
        values=np.array([statistics.median(r['times_ms'][stage] for r in profiles if r['implementation']==name)
                         for name in ['FrozenNeuralPoisson','CachedFrozenNeuralPoisson']])
        axes[1].bar([0,1],values,bottom=bottom,color=color,label=label);bottom+=values
    axes[1].set(xticks=[0,1],xticklabels=['Original','Cached coefficients'],ylabel='Median phase wall time (ms)',
                title='One 512-cell verification: cost breakdown')
    axes[1].legend(fontsize=7)
    for ax in axes:ax.grid(axis='y',alpha=.2)
    fig.tight_layout();fig.savefig(out/'cost_comparison.png',dpi=160);plt.close(fig)
    source_files=[Path(__file__),ROOT/'scripts/interval_cached.py',ROOT/'scripts/interval_certification.py',ROOT/'scripts/certification_core.py']
    metadata=dict(case_id='poisson_1d_sine_k1_seed7',case_sha256=sha256(CASE/'case.json'),
                  repeats=args.repeats,interval_precision_digits=30,cap=2048,initial_cells=32,
                  targets=[str(t) for t in TARGETS],python=sys.version,platform=platform.platform(),
                  numpy=np.__version__,scipy=scipy.__version__,mpmath=mpmath.__version__,
                  thread_settings={k:os.environ.get(k) for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']},
                  included=['case/weight loading and validation','coefficient caching when used','float pilot when used',
                            'all unsuccessful mesh trials','grid generation and residual/derivative evaluation',
                            'integration, continuous extrema bounds, stopping tests and safe endpoint export'],
                  excluded=['Python startup and module imports','warmup and audit checks',
                            'exact-solution reference evaluations','JSON/figure disk writes and console logging'],
                  timing_method='Sequential, mode order rotated per repeat, warmed implementations, fresh model per request; wall and process CPU time',
                  sources={str(p.relative_to(ROOT)):sha256(p) for p in source_files},
                  checks_passed=len(checks),limitations=['single small smooth 1D case','three repeats by default; preliminary measurements',
                                                      'Float64 is not an outward-rounded certificate',
                                                      'Interval guarantees retain the audited backend/math assumptions',
                                                      'FD is a scoped cost comparison, not a substitute guarantee for this PINN'])
    (out/'config.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print('SUMMARY',json.dumps(summary),flush=True)
    print('RESULT_DIRECTORY',out,flush=True)


if __name__=='__main__':main()
