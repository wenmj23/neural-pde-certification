"""Evaluate frozen seed49 networks and newly parameterized controlled peaks.

All timings include mesh selection and all certificate attempts; reference
checks and diagnostic-only probes are outside timings. No policy sees the
reference value or the case label. Outputs are written to a new directory.
"""
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path
import argparse, gc, hashlib, json, platform, statistics

import flint
from arb_mesh_selection import request as old_request, RecordedCubic, indicators
from arb_mesh_check_forecast import request as growth_request
from arb_cubic_reconstruction import SixthOrderDeep, BumpSix, SineSix, CubicVerifier
from arb_degree_selection import fraction, pack
from mesh_problem_families import GaussianMixture, MPGaussianReference

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT/'configs/frozen_networks_holdout_2026-10-30.json'
REFERENCE = None
POLICIES = ('cubic_uniform16', 'auto_heap16', 'nested_left125_16', 'nested_left_stride2_125_16', 'nested_left_stride3_125_16', 'nested_left_fcguard_125_16', 'nested_left_fcfast_125_16')
TARGETS = (F(51, 50), F(201, 200))

def main():
    out = ROOT/'results/generated/48_holdout_study'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True, exist_ok=False)
    fixture_data=json.loads(FIXTURE.read_text())
    gen_dir=ROOT/fixture_data['config']['generator_run']
    gen=json.loads((gen_dir/'config.json').read_text())
    if gen['status']!='complete':raise ValueError('New case generation incomplete')
    if hashlib.sha256(FIXTURE.read_bytes()).hexdigest()!=gen['fixture_sha256']:raise ValueError('New fixture mismatch')
    for name,expected in gen['source_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=expected:raise ValueError('Prelocked source changed: '+name)
    if list(POLICIES)!=gen['policies'] or [str(x) for x in TARGETS]!=gen['targets']:raise ValueError('Prelocked policy/target mismatch')
    refs_path=gen_dir/'network_references.json'
    sources=[ROOT/name for name in gen['source_sha256']]+[FIXTURE,gen_dir/'config.json',refs_path]
    snapshot = out/'source_snapshot'; snapshot.mkdir()
    for src in sources:
        if src.suffix == '.py': (snapshot/src.name).write_bytes(src.read_bytes())
    config = dict(status='running', suite='prelocked_seed49_and_new_shapes', cases=[], targets=[str(t) for t in TARGETS], policies=list(POLICIES), repeats=3,
                  max_intervals=4096, seconds=30, source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                  python=platform.python_version(), flint=flint.__version__,
                  timing='Sequential; fresh cache per request; all attempts and selection included; imports, garbage collection, references, diagnostic probes and file writes excluded.',
                  limitations='Four seed49 frozen networks generated after rules were locked; k2 is a new source frequency. Three new manufactured peak parameters were preselected. Same training and PDE family; not independent software audit. Forecast uses only completed certificate components; exact goal remains mandatory. Intermediate planned grids are not certified. Both targets are certificate enclosure ratios, not solution error.')
    def save(name, data):
        p=out/name; p.write_text(json.dumps(data,indent=2)+'\n'); assert json.loads(p.read_text()) == data
    def checkpoint(): save('config.json',config);save('runs.json',runs);save('checks.json',checks)
    fixture=FIXTURE; sha=hashlib.sha256(fixture.read_bytes()).hexdigest()
    cases=[(row['label'],lambda label=row['label']: SixthOrderDeep(fixture,label,sha)) for row in fixture_data['networks']]
    cases += [(label,lambda terms=terms:GaussianMixture(terms)) for label,terms in gen['shapes'].items()]
    refs=json.loads(refs_path.read_text())
    runs=[]; checks=[]; diagnostics=[]
    config['cases']=[x[0] for x in cases]; checkpoint()
    try:
        for label,terms in gen['shapes'].items():
            ref=MPGaussianReference(terms).maximum(relative=F(1,10**8),absolute=F(1,10**14),cap=5000)
            if not ref['converged']: raise ArithmeticError('New shape continuous reference failed '+label)
            refs[label]=ref
            save('shape_references.json',{k:refs[k] for k in gen['shapes'] if k in refs})
            print('SHAPE_REFERENCE',label,'complete',flush=True)
        for label,factory in cases:
            v=RecordedCubic(factory(),'moment');v.uniform(16);result=v.certify();_,c=indicators(v)
            diagnostics.append({'label':label,'effective_cell_fraction':c,'initial_enclosure':pack(result)})
        save('diagnostics.json',diagnostics)
        for repeat in range(3):
            for label,factory in cases:
                for target in TARGETS:
                    shift=(2*repeat)%len(POLICIES); order=POLICIES[shift:]+POLICIES[:shift]
                    for policy in order:
                        gc.collect(); r=growth_request(factory,policy,target,cap=4096,seconds=30)
                        r.update(label=label,repeat=repeat);runs.append(r)
                        if label in refs: lo,hi=fraction(refs[label]['lower']),fraction(refs[label]['upper'])
                        else: raise AssertionError('Missing independent reference '+label)
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
            if r['policy'].startswith('cubic_uniform') or not r['final']: continue
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
