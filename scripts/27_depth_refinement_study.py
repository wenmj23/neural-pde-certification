"""Follow up study25 with cheaper centering and derivative-only subdivision."""
from datetime import datetime,timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import statistics
from arb_deep_network import DeepPoisson
from arb_deep_refinements import RefinedDeepPoisson
from arb_shaped_remainder import ShapedVerifier

ROOT=Path(__file__).resolve().parents[1]
DRIVER=ROOT/'scripts/25_depth_certification_study.py'
spec=importlib.util.spec_from_file_location('study25',DRIVER)
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
FIXTURE=base.FIXTURE


def factory(path,label,sha,method):
    cls=DeepPoisson if method in ('chain','centered') else RefinedDeepPoisson
    return cls(path,label,sha,method)


def main():
    out=ROOT/'results/generated/27_depth_refinements'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True,exist_ok=False)
    ref_run=max(p.parent for p in (ROOT/'results/generated/24_depth_reference').glob('*/config.json') if json.loads(p.read_text()).get('status')=='complete')
    refs=json.loads((ref_run/'references.json').read_text());rows=json.loads(FIXTURE.read_text())['networks'];sha=hashlib.sha256(FIXTURE.read_bytes()).hexdigest();checks=[]
    methods=['chain','centered','value_centered','centered_lite','split']
    sources=[Path(__file__),DRIVER,ROOT/'scripts/arb_deep_network.py',ROOT/'scripts/arb_deep_refinements.py',ROOT/'scripts/arb_shaped_remainder.py',ROOT/'scripts/arb_reconstruction.py',FIXTURE,ref_run/'references.json']
    (out/'source_snapshot').mkdir()
    for p in sources:
        if p.suffix=='.py':(out/'source_snapshot'/p.name).write_bytes(p.read_bytes())
    config=dict(status='running',repeats=3,methods=methods,target='51/50',reconstruction='quadratic_moment',
                sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                scope='Same 18 frozen networks; derivative-only split adds interval calls but no residual interpolation nodes; all costs included',
                timing='Reuse study25 complete-request driver with a factory dispatch; 8s soft cap checked between full mesh attempts, cap N=2048; rotated method order.')
    def save(name,value):(out/name).write_text(json.dumps(value,indent=2)+'\n')
    save('config.json',config)
    def check(name,ok):
        checks.append(dict(name=name,passed=bool(ok)));save('checks.json',checks)
        if not ok:raise AssertionError(name)
    def consistent(label,s,name):check(name,base.frac(s['lower'])<=base.frac(refs[label]['lower']) and base.frac(s['upper'])>=base.frac(refs[label]['upper']))
    fixed=[]
    for row in rows:
        label=row['label']
        for n in (64,256):
            for method in methods:
                p=factory(FIXTURE,label,sha,method);v=ShapedVerifier(p,'moment');v.uniform(n);s=base.serialize(v.certify())
                consistent(label,s,f'fixed_ref_{label}_{n}_{method}')
                check(f'sample_budget_{label}_{n}_{method}',s['residual_evaluations']==2*n+1)
                fixed.append(dict(label=label,derivatives=method,**s,network_forward_calls=p.forward_calls,
                                  extra_center_passes=getattr(p,'center_passes',None),derivative_bound_subcells=getattr(p,'bound_subcells',0)))
        save('fixed_grid.json',fixed);print('FOLLOWUP_FIXED',label,flush=True)
    # The imported request function resolves this module variable at runtime;
    # source files are unchanged and snapshots retain the executed versions.
    base.DeepPoisson=factory;runs=[]
    for repeat in range(3):
        for row in rows:
            ordered=methods[repeat:]+methods[:repeat]
            for method in ordered:
                r=base.request(row['label'],sha,method,'quadratic_moment');r.update(repeat=repeat,depth=row['depth'],frequency=row['frequency'])
                if r['final'] is not None:consistent(row['label'],r['final'],f'timed_ref_{repeat}_{row["label"]}_{method}')
                runs.append(r);save('runs.json',runs)
                print('FOLLOWUP_TIMED',repeat,row['label'],method,r['status'],'N',r['final']['intervals'] if r['final'] else None,'ms',round(r['wall_ms'],2),flush=True)
    summary=[]
    for row in rows:
        for method in methods:
            group=[r for r in runs if r['label']==row['label'] and r['derivatives']==method]
            summary.append(dict(label=row['label'],depth=row['depth'],frequency=row['frequency'],derivatives=method,successes=sum(r['success'] for r in group),
                                intervals=group[0]['final']['intervals'] if group[0]['final'] else None,
                                median_wall_ms=statistics.median(r['wall_ms'] for r in group),min_wall_ms=min(r['wall_ms'] for r in group),max_wall_ms=max(r['wall_ms'] for r in group)))
    save('summary.json',summary);config.update(status='complete',checks_passed=len(checks),fixed_configurations=len(fixed),timed_requests=len(runs),successful_requests=sum(r['success'] for r in runs));save('config.json',config)
    print('COMPLETED',out,'checks',len(checks),'success',config['successful_requests'],'/',len(runs),flush=True)


if __name__=='__main__':main()
