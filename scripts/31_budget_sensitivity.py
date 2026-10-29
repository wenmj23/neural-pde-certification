"""Check cubic/quad timing conclusions under a finer common budget schedule."""
from datetime import datetime,timezone
from fractions import Fraction as F
from math import ceil
import gc
import hashlib
import importlib.util
import json
from pathlib import Path
import statistics
import time
from arb_cubic_reconstruction import SixthOrderDeep,CubicVerifier
from arb_shaped_remainder import ShapedVerifier

ROOT=Path(__file__).resolve().parents[1];FIXTURE=ROOT/'configs/frozen_networks_depth_2026-10-29.json'
DRIVER=ROOT/'scripts/25_depth_certification_study.py'
spec=importlib.util.spec_from_file_location('study25',DRIVER);base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)


def request(label,sha,degree):
    gc.collect();start=time.perf_counter();p=SixthOrderDeep(FIXTURE,label,sha)
    v=ShapedVerifier(p,'moment') if degree==2 else CubicVerifier(p,'moment')
    m=4;history=[];last=None;success=False;status='mesh_cap'
    while m<=1024:
        v.uniform(6*m//degree);last=base.serialize(v.certify());lo,hi=base.frac(last['lower']),base.frac(last['upper']);success=lo>0 and hi<=F(51,50)*lo
        history.append(dict(m=m,sample_budget=6*m+1,lower=last['lower'],upper=last['upper']))
        if success:status='success';break
        if time.perf_counter()-start>20:status='time_cap';break
        if m==1024:break
        m=min(1024,ceil(5*m/4))
    return dict(label=label,degree=degree,status=status,success=bool(success),final=last,history=history,wall_ms=(time.perf_counter()-start)*1000)


def main():
    out=ROOT/'results/generated/31_budget_sensitivity'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');out.mkdir(parents=True,exist_ok=False)
    sha=hashlib.sha256(FIXTURE.read_bytes()).hexdigest();labels=['depth1_k1_seed7','depth2_k1_seed27']+[f'depth3_k{k}_seed{s}' for k in (1,5) for s in (7,17,27)]
    ref_run=max(p.parent for p in (ROOT/'results/generated/24_depth_reference').glob('*/config.json') if json.loads(p.read_text()).get('status')=='complete');refs=json.loads((ref_run/'references.json').read_text())
    sources=[Path(__file__),DRIVER,ROOT/'scripts/arb_cubic_reconstruction.py',ROOT/'scripts/arb_deep_network.py',ROOT/'scripts/arb_shaped_remainder.py',ROOT/'scripts/arb_reconstruction.py',FIXTURE,ref_run/'references.json']
    (out/'source_snapshot').mkdir()
    for p in sources:
        if p.suffix=='.py':(out/'source_snapshot'/p.name).write_bytes(p.read_bytes())
    config=dict(status='running',repeats=3,labels=labels,target='51/50',growth='m <- ceil(5m/4), common B=6m+1, start m=4',
                timing='All fresh model, nonnested budget attempts and exports included; references and file writes excluded; 20s soft cap checked between meshes.',
                scope='Selected sensitivity subset, not held-out independent evidence; no use of known solution for budget choice',
                sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    def save(name,value):(out/name).write_text(json.dumps(value,indent=2)+'\n')
    save('config.json',config);checks=[];runs=[]
    for repeat in range(3):
        for label in labels:
            for degree in ([2,3] if repeat%2==0 else [3,2]):
                r=request(label,sha,degree);r['repeat']=repeat;ref=refs[label]
                ok=base.frac(r['final']['lower'])<=base.frac(ref['lower']) and base.frac(r['final']['upper'])>=base.frac(ref['upper'])
                checks.append(dict(name=f'reference_{repeat}_{label}_{degree}',passed=bool(ok)));save('checks.json',checks)
                if not ok:raise AssertionError('Reference not contained')
                runs.append(r);save('runs.json',runs)
                print('FINE_BUDGET',repeat,label,degree,r['status'],'N',r['final']['intervals'],'ms',round(r['wall_ms'],1),flush=True)
    summary=[]
    for label in labels:
        for degree in (2,3):
            group=[r for r in runs if r['label']==label and r['degree']==degree]
            summary.append(dict(label=label,degree=degree,successes=sum(r['success'] for r in group),intervals=group[0]['final']['intervals'],
                                median_wall_ms=statistics.median(r['wall_ms'] for r in group),min_wall_ms=min(r['wall_ms'] for r in group),max_wall_ms=max(r['wall_ms'] for r in group)))
    save('summary.json',summary);config.update(status='complete',checks_passed=len(checks),timed_requests=len(runs),successful_requests=sum(r['success'] for r in runs));save('config.json',config)
    print('COMPLETED',out,'success',config['successful_requests'],'/',len(runs),flush=True)


if __name__=='__main__':main()
