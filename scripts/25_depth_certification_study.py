"""Depth benchmark: equal residual samples and fresh matched-gap requests."""
from datetime import datetime,timezone
from fractions import Fraction
import gc
import hashlib
import json
from pathlib import Path
import statistics
import time

from flint import arb
from arb_deep_network import DeepPoisson
from arb_reconstruction import Verifier,norm_enclosure
from arb_shaped_remainder import ShapedVerifier
from arb_certification import exact_fraction,export

ROOT=Path(__file__).resolve().parents[1]
FIXTURE=ROOT/'configs/frozen_networks_depth_2026-10-29.json'


def frac(e):return Fraction(int(e['numerator']),int(e['denominator']))


def serialize(r):
    fields=['lower','upper','z_lower','z_upper','remainder_upper']
    if 'constant_remainder_upper' in r:fields.append('constant_remainder_upper')
    return dict(**{k:export(r[k],'down' if k in ('lower','z_lower') else 'up') for k in fields},
                intervals=r['intervals'],residual_evaluations=r['residual_evaluations'],
                derivative_interval_evaluations=r['curvature_evaluations'],degree=r['degree'])


def global_residual_upper(cells):
    maximum=arb(0)
    for c in cells:
        coeff=[c['r0'],2*c['rm']-(c['r0']+c['r1'])/2,c['r1']]
        _,hi=norm_enclosure(coeff,2)
        maximum=max(maximum,(hi+c['delta']).upper())
    return (maximum/8).upper()


def make_verifier(p,reconstruction):
    if reconstruction=='linear':return Verifier(p,degree=1,subdivision=2)
    return ShapedVerifier(p,mode='constant' if reconstruction=='quadratic_constant' else 'moment')


def request(label,sha,derivatives,reconstruction,target=Fraction(51,50),cap=2048,seconds=8):
    gc.collect();start=time.perf_counter();cpu=time.process_time()
    p=DeepPoisson(FIXTURE,label,sha,derivatives);v=make_verifier(p,reconstruction)
    budget=16;history=[];last=None;success=False;status='mesh_cap'
    while budget<=cap:
        n=2*budget if reconstruction=='linear' else budget
        v.uniform(n)
        try:
            last=serialize(v.certify());lo,hi=frac(last['lower']),frac(last['upper'])
            success=lo>0 and hi<=target*lo
            history.append(dict(budget=budget,intervals=n,status='finite',lower=last['lower'],upper=last['upper']))
        except ArithmeticError as error:
            history.append(dict(budget=budget,intervals=n,status='safe_rejection',reason=str(error)))
        if success:status='success';break
        if time.perf_counter()-start>=seconds:status='time_cap';break
        budget*=2
    wall=(time.perf_counter()-start)*1000;cpu=(time.process_time()-cpu)*1000
    return dict(label=label,derivatives=derivatives,reconstruction=reconstruction,target=str(target),
                status=status,success=success,final=last,history=history,wall_ms=wall,cpu_ms=cpu,
                network_forward_calls=p.forward_calls,
                note='Time cap checked between complete mesh attempts; final attempt can exceed soft time cap.')


def main():
    out=ROOT/'results/generated/25_depth_certification'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True,exist_ok=False)
    reference_run=max(p.parent for p in (ROOT/'results/generated/24_depth_reference').glob('*/config.json')
                      if json.loads(p.read_text()).get('status')=='complete')
    references=json.loads((reference_run/'references.json').read_text());data=json.loads(FIXTURE.read_text());sha=hashlib.sha256(FIXTURE.read_bytes()).hexdigest()
    sources=[Path(__file__),ROOT/'scripts/arb_deep_network.py',ROOT/'scripts/arb_reconstruction.py',ROOT/'scripts/arb_shaped_remainder.py',ROOT/'scripts/arb_certification.py',ROOT/'scripts/arb_problem_families.py',FIXTURE,reference_run/'references.json']
    (out/'source_snapshot').mkdir()
    for p in sources:
        if p.suffix=='.py':(out/'source_snapshot'/p.name).write_bytes(p.read_bytes())
    config=dict(status='running',sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                reference=str(reference_run.relative_to(ROOT)),bits=103,repeats=3,target='51/50',
                soft_request_time_cap_s=8,quadratic_mesh_cap=2048,linear_mesh_cap=4096,
                fixed_grid_protocol='Each quadratic N matches linear 2N: both 2N+1 residual samples; derivative counts and postprocessing differ.',
                timing='Fresh model/cache; rotated methods; all failed meshes, model loading, integration, norm bounds and endpoint export included; reference checks/file writes excluded.',
                scope='Unit-interval zero-Dirichlet Poisson with frozen 1/2/3-hidden-layer tanh candidates; conditional interval proof, no independent/formal software audit')
    def save(name,value):(out/name).write_text(json.dumps(value,indent=2)+'\n')
    save('config.json',config);checks=[]
    def check(name,ok):
        checks.append(dict(name=name,passed=bool(ok)));save('checks.json',checks)
        if not ok:raise AssertionError(name)
    def consistent(label,s,name):
        ref=references[label]
        check(name,frac(s['lower'])<=frac(ref['lower']) and frac(s['upper'])>=frac(ref['upper']))
    fixed=[]
    for row in data['networks']:
        label=row['label']
        for derivatives in ('exp','chain','centered'):
            for n in (16,64,256):
                p=DeepPoisson(FIXTURE,label,sha,derivatives);v=ShapedVerifier(p,'constant');v.uniform(n)
                try:
                    constant=serialize(v.certify());counts=(v.residual_evaluations,v.curvature_evaluations)
                    residual_max=export(global_residual_upper(v.prepare(v.nodes)),'up')
                    v.mode='moment';moment=serialize(v.certify())
                    check(f'reused_data_{label}_{derivatives}_{n}',counts==(v.residual_evaluations,v.curvature_evaluations))
                    for method,s in [('quadratic_constant',constant),('quadratic_moment',moment)]:
                        check(f'count_{label}_{derivatives}_{n}_{method}',s['residual_evaluations']==2*n+1 and s['derivative_interval_evaluations']==n)
                        consistent(label,s,f'fixed_ref_{label}_{derivatives}_{n}_{method}')
                        fixed.append(dict(label=label,derivatives=derivatives,reconstruction=method,status='finite',sample_budget=2*n+1,**s,global_residual_upper=residual_max))
                    check(f'moment_never_wider_{label}_{derivatives}_{n}',frac(moment['upper'])<=frac(constant['upper']) and frac(moment['lower'])>=frac(constant['lower']))
                    check(f'max_residual_reference_{label}_{derivatives}_{n}',frac(residual_max)>=frac(references[label]['upper']))
                except ArithmeticError as error:
                    fixed.append(dict(label=label,derivatives=derivatives,reconstruction='quadratic_pair',status='safe_rejection',sample_budget=2*n+1,reason=str(error)))
                p=DeepPoisson(FIXTURE,label,sha,derivatives);v=Verifier(p,1,2);v.uniform(2*n)
                try:
                    s=serialize(v.certify());consistent(label,s,f'linear_ref_{label}_{derivatives}_{n}')
                    check(f'linear_count_{label}_{derivatives}_{n}',s['residual_evaluations']==2*n+1 and s['derivative_interval_evaluations']==2*n)
                    fixed.append(dict(label=label,derivatives=derivatives,reconstruction='linear',status='finite',sample_budget=2*n+1,**s))
                except ArithmeticError as error:
                    fixed.append(dict(label=label,derivatives=derivatives,reconstruction='linear',status='safe_rejection',sample_budget=2*n+1,reason=str(error)))
        save('fixed_grid.json',fixed)
        print('SCREENED',label,flush=True)
    methods=[('exp','quadratic_moment'),('chain','linear'),('chain','quadratic_constant'),('chain','quadratic_moment'),
             ('centered','linear'),('centered','quadratic_constant'),('centered','quadratic_moment')]
    runs=[]
    for repeat in range(3):
        for row in data['networks']:
            label=row['label'];ordered=methods[repeat:]+methods[:repeat]
            for derivatives,reconstruction in ordered:
                r=request(label,sha,derivatives,reconstruction);r.update(repeat=repeat,depth=row['depth'],frequency=row['frequency'],seed=row['seed'])
                if r['final'] is not None:consistent(label,r['final'],f'timed_ref_{repeat}_{label}_{derivatives}_{reconstruction}')
                if r['success']:check(f'target_{repeat}_{label}_{derivatives}_{reconstruction}',frac(r['final']['upper'])<=Fraction(51,50)*frac(r['final']['lower']))
                runs.append(r);save('runs.json',runs)
                print('TIMED',repeat,label,derivatives,reconstruction,r['status'],
                      'N',r['final']['intervals'] if r['final'] else None,'ms',round(r['wall_ms'],1),flush=True)
    summary=[]
    for row in data['networks']:
        for derivatives,reconstruction in methods:
            group=[r for r in runs if (r['label'],r['derivatives'],r['reconstruction'])==(row['label'],derivatives,reconstruction)]
            summary.append(dict(label=row['label'],depth=row['depth'],frequency=row['frequency'],seed=row['seed'],derivatives=derivatives,reconstruction=reconstruction,
                                successes=sum(r['success'] for r in group),intervals=group[0]['final']['intervals'] if group[0]['final'] else None,
                                median_wall_ms=statistics.median(r['wall_ms'] for r in group),min_wall_ms=min(r['wall_ms'] for r in group),max_wall_ms=max(r['wall_ms'] for r in group)))
    save('summary.json',summary)
    config.update(status='complete',checks_passed=len(checks),fixed_records=len(fixed),fixed_safe_rejections=sum(r['status']=='safe_rejection' for r in fixed),timed_requests=len(runs),successful_requests=sum(r['success'] for r in runs))
    save('config.json',config)
    print('COMPLETED',out,'checks',len(checks),'successful',config['successful_requests'],'/',len(runs),flush=True)


if __name__=='__main__':main()
