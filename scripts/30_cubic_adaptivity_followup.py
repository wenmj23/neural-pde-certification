"""Combine higher-order reconstruction with existing local refinement."""
from datetime import datetime,timezone
from fractions import Fraction as F
import gc
import hashlib
import importlib.util
import json
from pathlib import Path
import statistics
import time
from arb_cubic_reconstruction import CubicVerifier,SixthOrderDeep,SineSix,BumpSix
from arb_shaped_remainder import ShapedVerifier

ROOT=Path(__file__).resolve().parents[1]
FIXTURE=ROOT/'configs/frozen_networks_depth_2026-10-29.json'
DRIVER=ROOT/'scripts/25_depth_certification_study.py'
spec=importlib.util.spec_from_file_location('study25',DRIVER);base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)


def request(factory,degree,adaptive,target):
    gc.collect();start=time.perf_counter();p=factory()
    v=ShapedVerifier(p,'moment') if degree==2 else CubicVerifier(p,'moment')
    m=8;v.uniform(6*m//degree);history=[];last=None;status='mesh_cap';success=False
    while m<=1024:
        n=6*m//degree
        if adaptive:v.refine_to(n)
        else:v.uniform(n)
        try:
            last=base.serialize(v.certify());lo,hi=base.frac(last['lower']),base.frac(last['upper']);success=lo>0 and hi<=target*lo
            history.append(dict(sample_budget=6*m+1,status='finite',lower=last['lower'],upper=last['upper']))
        except ArithmeticError as error:history.append(dict(sample_budget=6*m+1,status='safe_rejection',reason=str(error)))
        if success:status='success';break
        if time.perf_counter()-start>8:status='time_cap';break
        m*=2
    wall=(time.perf_counter()-start)*1000
    return dict(degree=degree,adaptive=adaptive,target=str(target),status=status,success=bool(success),final=last,history=history,wall_ms=wall)


def main():
    out=ROOT/'results/generated/30_cubic_adaptivity'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');out.mkdir(parents=True,exist_ok=False)
    sha=hashlib.sha256(FIXTURE.read_bytes()).hexdigest();refs_run=max(p.parent for p in (ROOT/'results/generated/24_depth_reference').glob('*/config.json') if json.loads(p.read_text()).get('status')=='complete');refs=json.loads((refs_run/'references.json').read_text())
    factories={'sine50':lambda:SineSix(50),'bump':lambda:BumpSix('0.005')}
    for label in ['depth3_k1_seed7','depth3_k5_seed7']:factories[label]=lambda label=label:SixthOrderDeep(FIXTURE,label,sha)
    sources=[Path(__file__),DRIVER,ROOT/'scripts/arb_cubic_reconstruction.py',ROOT/'scripts/arb_deep_network.py',ROOT/'scripts/arb_shaped_remainder.py',ROOT/'scripts/arb_reconstruction.py',FIXTURE,refs_run/'references.json']
    (out/'source_snapshot').mkdir()
    for p in sources:
        if p.suffix=='.py':(out/'source_snapshot'/p.name).write_bytes(p.read_bytes())
    config=dict(status='running',repeats=3,protocol='Start B=49, then B-1 doubles; quadratic N=(B-1)/2, cubic N=(B-1)/3; same initial uniform mesh within each degree before optional refinement; all search costs included.',
                scope='Four selected follow-up tasks, not a new independent dataset; interval assumptions unchanged',
                sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    def save(name,value):(out/name).write_text(json.dumps(value,indent=2)+'\n')
    save('config.json',config);checks=[];runs=[]
    methods=[(2,False),(2,True),(3,False),(3,True)]
    tasks=[('sine50',F(6,5)),('bump',F(51,50)),('depth3_k1_seed7',F(51,50)),('depth3_k5_seed7',F(51,50))]
    for repeat in range(3):
        for label,target in tasks:
            for degree,adaptive in methods[repeat:]+methods[:repeat]:
                r=request(factories[label],degree,adaptive,target);r.update(label=label,repeat=repeat)
                if r['final'] is not None:
                    lo,hi=(base.frac(refs[label]['lower']),base.frac(refs[label]['upper'])) if label in refs else (F(1,10000),F(1,10000))
                    ok=base.frac(r['final']['lower'])<=lo and base.frac(r['final']['upper'])>=hi
                    checks.append(dict(name=f'reference_{repeat}_{label}_{degree}_{adaptive}',passed=bool(ok)));save('checks.json',checks)
                    if not ok:raise AssertionError('Reference not contained')
                runs.append(r);save('runs.json',runs)
                print('ADAPT_CUBIC',repeat,label,degree,adaptive,r['status'],'N',r['final']['intervals'] if r['final'] else None,'ms',round(r['wall_ms'],2),flush=True)
    summary=[]
    for label,target in tasks:
        for degree,adaptive in methods:
            group=[r for r in runs if r['label']==label and r['degree']==degree and r['adaptive']==adaptive]
            summary.append(dict(label=label,target=str(target),degree=degree,adaptive=adaptive,successes=sum(r['success'] for r in group),
                                intervals=group[0]['final']['intervals'] if group[0]['final'] else None,
                                median_wall_ms=statistics.median(r['wall_ms'] for r in group),min_wall_ms=min(r['wall_ms'] for r in group),max_wall_ms=max(r['wall_ms'] for r in group)))
    save('summary.json',summary);config.update(status='complete',checks_passed=len(checks),timed_requests=len(runs),successful_requests=sum(r['success'] for r in runs));save('config.json',config)
    print('COMPLETED',out,'success',config['successful_requests'],'/',len(runs),flush=True)


if __name__=='__main__':main()
