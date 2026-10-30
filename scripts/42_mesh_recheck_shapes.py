"""Experimental 32-cell mesh recheck with frozen networks and controls.

All timings include mesh selection and all certificate attempts; reference
checks and diagnostic-only probes are outside timings. No policy sees the
reference value or the case label. Outputs are written to a new directory.
"""
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path
import gc, hashlib, json, platform, statistics

import flint
from arb_mesh_selection import request as old_request, RecordedCubic, indicators
from arb_mesh_recheck import request as recheck_request
from arb_cubic_reconstruction import SixthOrderDeep, BumpSix, SineSix, CubicVerifier
from arb_degree_selection import fraction, pack
from mesh_problem_families import GaussianMixture

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT/'configs/frozen_networks_depth_2026-10-29.json'
REFERENCE = ROOT/'results/generated/38_mesh_selection/20261029T115506215623Z/references.json'
POLICIES = ('cubic_uniform16', 'auto_heap16', 'auto_heap_recheck16')
TARGETS = (F(51, 50), F(201, 200))

def main():
    out = ROOT/'results/generated/42_mesh_recheck_shapes'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True, exist_ok=False)
    sources = [Path(__file__), ROOT/'scripts/arb_mesh_selection.py', ROOT/'scripts/arb_mesh_recheck.py', ROOT/'scripts/mesh_problem_families.py', ROOT/'scripts/arb_cubic_reconstruction.py', ROOT/'scripts/arb_degree_selection.py', FIXTURE, REFERENCE]
    config = dict(status='running', cases=[], targets=[str(t) for t in TARGETS], policies=list(POLICIES), repeats=3,
                  max_intervals=4096, seconds=30, source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                  python=platform.python_version(), flint=flint.__version__,
                  timing='Sequential; fresh cache per request; all attempts and selection included; imports, garbage collection, references, diagnostic probes and file writes excluded.',
                  limitations='Five previously studied manufactured off-center, near-boundary and multipeak cases; no networks; no formal audit. Both targets are certificate enclosure ratios, not solution error.')
    def save(name, data):
        p=out/name; p.write_text(json.dumps(data,indent=2)+'\n'); assert json.loads(p.read_text()) == data
    def checkpoint(): save('config.json',config);save('runs.json',runs);save('checks.json',checks)
    sha=hashlib.sha256(FIXTURE.read_bytes()).hexdigest()
    shapes = {
        'shift_0137': [dict(center='0.137',width='0.007',weight='1')],
        'shift_0733': [dict(center='0.733',width='0.013',weight='1')],
        'boundary_097': [dict(center='0.97',width='0.004',weight='1')],
        'double_012_088': [dict(center='0.12',width='0.006',weight='1/2'),dict(center='0.88',width='0.009',weight='1/2')],
        'triple_015_047_083': [dict(center='0.15',width='0.008',weight='1/3'),dict(center='0.47',width='0.011',weight='1/3'),dict(center='0.83',width='0.006',weight='1/3')],
    }
    cases=[(label,lambda terms=terms:GaussianMixture(terms)) for label,terms in shapes.items()]
    config['shapes']=shapes
    refs=json.loads(REFERENCE.read_text()); runs=[]; checks=[]; diagnostics=[]
    config['cases']=[x[0] for x in cases]; checkpoint()
    try:
        for label,factory in cases:
            v=RecordedCubic(factory(),'moment');v.uniform(16);result=v.certify();_,c=indicators(v)
            diagnostics.append({'label':label,'effective_cell_fraction':c,'initial_enclosure':pack(result)})
        save('diagnostics.json',diagnostics)
        for repeat in range(3):
            for label,factory in cases:
                for target in TARGETS:
                    order=POLICIES[repeat:]+POLICIES[:repeat]
                    for policy in order:
                        gc.collect(); fn=recheck_request if policy=='auto_heap_recheck16' else old_request; r=fn(factory,policy,target,cap=4096,seconds=30)
                        r.update(label=label,repeat=repeat);runs.append(r)
                        if label in refs: lo,hi=fraction(refs[label]['lower']),fraction(refs[label]['upper'])
                        else: raise AssertionError('shifted shape requires independent reference')
                        if r['final']:
                            final=r['final'];check=fraction(final['lower'])<=lo<=hi<=fraction(final['upper'])
                            checks.append({'label':label,'target':str(target),'policy':policy,'repeat':repeat,'reference_contained':check})
                            if not check: raise AssertionError('reference not contained')
                        checkpoint()
                        print(label,str(target),policy,round(r['wall_ms'],2),r['success'],flush=True)
        summary=[]
        for label,_ in cases:
            for target in TARGETS:
                for policy in POLICIES:
                    rows=[r for r in runs if r['label']==label and r['target']==str(target) and r['policy']==policy]
                    summary.append(dict(label=label,target=str(target),policy=policy,successes=sum(r['success'] for r in rows),
                                        median_wall_ms=statistics.median(r['wall_ms'] for r in rows),
                                        wall_ms=[r['wall_ms'] for r in rows],intervals=[r['final']['intervals'] if r['final'] else None for r in rows],
                                        routes=[r['locked_mesh'] for r in rows],recheck_switched=[any(t.get('recheck_switched_to_uniform',False) for t in r['trace']) for r in rows],work=[r['work'] for r in rows]))
        save('summary.json',summary)
        # Original-kernel exact equality on every distinct locally refined final mesh.
        matches=[]; seen=set()
        for r in runs:
            if r['policy']=='cubic_uniform16' or not r['final']: continue
            key=(r['label'],r['degree'],tuple(r['final_nodes']))
            if key in seen: continue
            seen.add(key);factory=dict(cases)[r['label']]
            original=pack(CubicVerifier(factory(),'moment').certify([F(x) for x in r['final_nodes']]))
            ok=original==r['final']; matches.append({'label':r['label'],'intervals':r['final']['intervals'],'matched':ok})
            if not ok: raise AssertionError('original-kernel mismatch')
        save('original_kernel_matches.json',matches)
        config.update(status='complete',requests=len(runs),successes=sum(r['success'] for r in runs),reference_checks=len(checks),original_meshes_checked=len(matches))
        checkpoint();print('COMPLETED',out,'requests',len(runs),'checks',len(checks),'original meshes',len(matches),flush=True)
    except Exception as exc:
        config.update(status='failed',error=repr(exc));checkpoint();raise

if __name__=='__main__': main()
