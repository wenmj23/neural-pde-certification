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

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT/'configs/frozen_networks_depth_2026-10-29.json'
REFERENCE = ROOT/'results/generated/24_depth_reference/20261029T044806546934Z/references.json'
POLICIES = ('cubic_uniform16', 'auto_heap16', 'auto_heap_recheck16')
TARGETS = (F(51, 50), F(201, 200))

def main():
    out = ROOT/'results/generated/41_mesh_recheck_study'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True, exist_ok=False)
    sources = [Path(__file__), ROOT/'scripts/arb_mesh_selection.py', ROOT/'scripts/arb_mesh_recheck.py', ROOT/'scripts/arb_cubic_reconstruction.py', ROOT/'scripts/arb_degree_selection.py', FIXTURE, REFERENCE]
    config = dict(status='running', cases=[], targets=[str(t) for t in TARGETS], policies=list(POLICIES), repeats=3,
                  max_intervals=4096, seconds=30, source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                  python=platform.python_version(), flint=flint.__version__,
                  timing='Sequential; fresh cache per request; all attempts and selection included; imports, garbage collection, references, diagnostic probes and file writes excluded.',
                  limitations='Twelve previously studied networks; seed27 not used in the immediate Oct30 diagnostic but seen in previous studies; two manufactured controls; no fresh-network holdout; no formal audit. Both targets are certificate enclosure ratios, not solution error.')
    def save(name, data):
        p=out/name; p.write_text(json.dumps(data,indent=2)+'\n'); assert json.loads(p.read_text()) == data
    def checkpoint(): save('config.json',config);save('runs.json',runs);save('checks.json',checks)
    sha=hashlib.sha256(FIXTURE.read_bytes()).hexdigest()
    cases=[]
    for depth in (1,3):
        for k in (1,5):
            for seed in (7,17,27):
                label=f'depth{depth}_k{k}_seed{seed}'
                cases.append((label, lambda label=label: SixthOrderDeep(FIXTURE,label,sha)))
    cases += [('narrow_peak',lambda:BumpSix('0.005')),('global_sine_50',lambda:SineSix(50))]
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
                        else: lo=hi=F(1,10000)
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
