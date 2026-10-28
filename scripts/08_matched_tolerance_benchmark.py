"""Select work using the certificate's own two-sided gap, not a reference error.

Ordinary float64 benchmark; include all refinement and failed-attempt costs.
"""
from pathlib import Path
import importlib.util,json,time,statistics
import numpy as np
from certification_core import SineError,BumpError,certificate,adaptive_nodes

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/generated/08_matched_tolerance_benchmark'


def load(name,file):
    spec=importlib.util.spec_from_file_location(name,Path(__file__).with_name(file))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def solve(problem,method,adaptive,ratio_target=1.2,cap=131072):
    start=time.perf_counter();intervals=32;nodes=None;queries=0;points=0;history=[]
    while intervals<=cap:
        if adaptive:
            nodes,meta=adaptive_nodes(problem,intervals,method,nodes)
            queries+=meta['bound_evaluations']
        else:nodes=np.linspace(0,1,intervals+1)
        result=certificate(problem,nodes,method);queries+=intervals;points+=len(nodes)
        gap=result['upper']/result['lower'] if result['lower']>0 else None
        success=(gap is not None and gap<=ratio_target) or result['upper']<=1e-12
        history.append(dict(intervals=intervals,lower=result['lower'],upper=result['upper'],gap=gap))
        if success:break
        intervals*=2
    return dict(success=bool(success),intervals=len(nodes)-1,lower=result['lower'],upper=result['upper'],
                gap=gap,total_ms=(time.perf_counter()-start)*1000,
                cumulative_residual_point_evaluations=points,
                cumulative_interval_bound_evaluations=queries,history=history)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    helper=load('methods07','07_high_frequency_methods.py')
    problem_list=[(SineError(20),1e-4),(SineError(50),1e-4),(BumpError(.005),1e-4)]
    problem_list+=list(helper.neural_problems())[:1]
    rows=[]
    for problem,reference in problem_list:
        for method,adaptive in [('constant',False),('linear',False),('linear',True)]:
            runs=[solve(problem,method,adaptive) for _ in range(5)]
            r=runs[0]
            assert all(v['success']==r['success'] and v['intervals']==r['intervals'] for v in runs)
            r.update(problem=problem.name,method=method,adaptive=adaptive,
                     reference_error=reference,ratio_to_reference=r['upper']/reference,
                     median_total_ms=statistics.median(v['total_ms'] for v in runs),
                     timing_samples_ms=[v['total_ms'] for v in runs])
            rows.append(r)
            print(json.dumps({k:r[k] for k in ['problem','method','adaptive','success','intervals','gap','median_total_ms','ratio_to_reference']}),flush=True)
    fd=load('checks06','06_known_error_and_cost_checks.py')
    fdrows=[]
    for n in (100,1000):
        ms=fd.median_ms(lambda:fd.finite_difference(n),30)
        x,u=fd.finite_difference(n);xe=np.linspace(0,1,10001)
        error=np.max(abs(np.interp(xe,x,u)-np.sin(np.pi*xe)))
        fdrows.append(dict(intervals=n,median_solve_ms=ms,reference_error=float(error)))
    (OUT/'records.json').write_text(json.dumps(rows,indent=2))
    (OUT/'finite_difference.json').write_text(json.dumps(fdrows,indent=2))
    print('FD scoped benchmark:',json.dumps(fdrows),flush=True)

    # Application-facing decisions use only upper/lower bounds and a prescribed
    # absolute tolerance. A failed acceptance test alone is not a rejection proof.
    decisions=[]
    for problem,ref in helper.neural_problems():
        for tolerance in (2e-4,5e-5):
            for method in ('constant','linear'):
                start=time.perf_counter();n=32;status='uncertain'
                while n<=16384:
                    result=certificate(problem,np.linspace(0,1,n+1),method)
                    if result['upper']<=tolerance:status='accept';break
                    if result['lower']>tolerance:status='reject';break
                    if n==16384:break
                    n*=2
                assert status!='accept' or ref<=tolerance+1e-12
                assert status!='reject' or ref>=tolerance-1e-10
                decisions.append(dict(problem=problem.name,tolerance=tolerance,method=method,
                    status=status,intervals=n,lower=result['lower'],upper=result['upper'],
                    reference_error=ref,total_ms=(time.perf_counter()-start)*1000))
    (OUT/'tolerance_decisions.json').write_text(json.dumps(decisions,indent=2))
    print('Tolerance decisions:',len(decisions),'reference-consistent decisions; statuses:',
          {s:sum(d['status']==s for d in decisions) for s in ['accept','reject','uncertain']},flush=True)


if __name__=='__main__':main()
