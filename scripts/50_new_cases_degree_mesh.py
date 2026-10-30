"""Matched quadratic/cubic and uniform/local study on frozen new cases.

This study is designed after seed49 generation; it is not a prelocked holdout.
All inputs, executions, budgets, certificates and failures are retained.
Reference values are only used after each timed request.
"""
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path
import gc, hashlib, json, platform, statistics
import flint

from arb_degree_selection import run_request, fraction, pack
from arb_mesh_selection import request as mesh_request
from arb_cubic_reconstruction import SixthOrderDeep, CubicVerifier
from arb_shaped_remainder import ShapedVerifier
from mesh_problem_families import GaussianMixture

ROOT=Path(__file__).resolve().parents[1]
FIXTURE=ROOT/'configs/frozen_networks_holdout_2026-10-30.json'
GEN=ROOT/'results/generated/47_new_cases/20261030T031028111446Z'
SHAPE_REFS=ROOT/'results/generated/48_holdout_study/20261030T031050542385Z/shape_references.json'
POLICIES=('fixed2','fixed3','quadratic_adaptive16','cubic_adaptive16','cubic_heap16','auto_heap16')
TARGETS=(F(51,50),F(201,200))

def main():
    fixture=json.loads(FIXTURE.read_text());g=json.loads((GEN/'config.json').read_text())
    if g['status']!='complete' or hashlib.sha256(FIXTURE.read_bytes()).hexdigest()!=g['fixture_sha256']:raise ValueError('Fixture mismatch')
    refs=json.loads((GEN/'network_references.json').read_text());refs.update(json.loads(SHAPE_REFS.read_text()))
    cases=[(row['label'],lambda label=row['label']:SixthOrderDeep(FIXTURE,label,g['fixture_sha256'])) for row in fixture['networks']]
    cases += [(label,lambda terms=terms:GaussianMixture(terms)) for label,terms in g['shapes'].items()]
    out=ROOT/'results/generated/50_new_cases_degree_mesh'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');out.mkdir(parents=True,exist_ok=False)
    names=['arb_degree_selection.py','arb_mesh_selection.py','arb_cubic_reconstruction.py','arb_shaped_remainder.py',
           'arb_reconstruction.py','arb_deep_network.py','arb_certification.py','arb_problem_families.py','mesh_problem_families.py',
           'interval_certification.py','18_continuous_reference_crosscheck.py']
    sources=[Path(__file__),*[ROOT/'scripts'/n for n in names],FIXTURE,GEN/'config.json',GEN/'network_references.json',SHAPE_REFS]
    snapshot=out/'source_snapshot';snapshot.mkdir()
    for p in sources:
        if p.suffix=='.py':(snapshot/p.name).write_bytes(p.read_bytes())
    config=dict(status='running',cases=[x[0] for x in cases],targets=[str(x) for x in TARGETS],policies=list(POLICIES),repeats=3,
                cap=4096,soft_seconds=30,max_attempts=40,flint=flint.__version__,python=platform.python_version(),
                source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                timing='Sequential; includes factory, caches, all attempts, mesh marking and exports; excludes imports, pre-request GC, output writes and validation-only references.',
                limitations='Four networks already generated and observed in preceding study; three manufactured shapes; no fresh holdout. Frozen candidate, zero BC, 1D Poisson and tanh; no formal audit.')
    runs=[];checks=[]
    def save(name,obj):
        p=out/name;p.write_text(json.dumps(obj,indent=2)+'\n');assert json.loads(p.read_text())==obj
    def checkpoint():save('config.json',config);save('runs.json',runs);save('checks.json',checks)
    checkpoint()
    try:
        for repeat in range(3):
            shift=2*repeat;order=POLICIES[shift:]+POLICIES[:shift]
            for label,factory in cases:
                ref=refs[label];lo,hi=fraction(ref['lower']),fraction(ref['upper'])
                for target in TARGETS:
                    for policy in order:
                        gc.collect()
                        if policy in ('fixed2','fixed3'):r=run_request(factory,policy,target,cap=4096,seconds=30,max_attempts=40)
                        else:r=mesh_request(factory,policy,target,cap=4096,seconds=30,max_attempts=40)
                        r.update(label=label,repeat=repeat);runs.append(r)
                        if r['final']:
                            a,b=fraction(r['final']['lower']),fraction(r['final']['upper']);passed=a<=lo<=hi<=b
                            checks.append({'label':label,'target':str(target),'policy':policy,'repeat':repeat,'reference_contained':passed})
                            if not passed:raise AssertionError('Reference not contained')
                        checkpoint();print(label,str(target),policy,round(r['wall_ms'],2),r['success'],flush=True)
        summary=[]
        for label,_ in cases:
            for target in TARGETS:
                for policy in POLICIES:
                    group=[r for r in runs if r['label']==label and r['target']==str(target) and r['policy']==policy]
                    summary.append(dict(label=label,target=str(target),policy=policy,successes=sum(x['success'] for x in group),
                                        median_ms=statistics.median(x['wall_ms'] for x in group),wall_ms=[x['wall_ms'] for x in group],
                                        final_intervals=[x['final']['intervals'] if x['final'] else None for x in group],
                                        work=[x['work'] for x in group]))
        save('summary.json',summary)
        seen={};matches=[]
        for r in runs:
            if r['final'] is None:continue
            if r['policy'] in ('fixed2','fixed3'):
                nodes=[F(i,r['final']['intervals']) for i in range(r['final']['intervals']+1)]
            else:nodes=[F(x) for x in r['final_nodes']]
            degree=r['final']['degree'];key=(r['label'],degree,tuple(nodes))
            if key not in seen:
                oracle=dict(cases)[r['label']]();v=ShapedVerifier(oracle,'moment') if degree==2 else CubicVerifier(oracle,'moment')
                seen[key]=pack(v.certify(nodes));matches.append({'label':r['label'],'degree':degree,'intervals':len(nodes)-1,'matched':seen[key]==r['final']})
            if seen[key]!=r['final']:raise AssertionError('Original kernel mismatch')
        save('original_kernel_matches.json',matches)
        config.update(status='complete',requests=len(runs),successful=sum(x['success'] for x in runs),reference_checks=len(checks),original_meshes_checked=len(matches));checkpoint()
    except Exception as exc:
        config.update(status='failed',error=repr(exc));checkpoint();raise
    print('COMPLETED',out,'requests',len(runs),'successes',config['successful'],'meshes',len(matches),flush=True)

if __name__=='__main__':main()
