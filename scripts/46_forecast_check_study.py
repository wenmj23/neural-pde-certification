"""Compare two forecast-based check schedules with fixed schedules.

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
from mesh_problem_families import GaussianMixture

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT/'configs/frozen_networks_depth_2026-10-29.json'
REFERENCE = ROOT/'results/generated/38_mesh_selection/20261029T115506215623Z/references.json'
POLICIES = ('cubic_uniform16', 'auto_heap16', 'nested_left125_16', 'nested_left_stride2_125_16', 'nested_left_stride3_125_16', 'nested_left_fcguard_125_16', 'nested_left_fcfast_125_16')
TARGETS = (F(51, 50), F(201, 200))

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--suite',choices=('development','evaluation'),required=True); args=parser.parse_args()
    out = ROOT/'results/generated/46_forecast_check_study'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True, exist_ok=False)
    sources = [Path(__file__), ROOT/'scripts/arb_mesh_selection.py', ROOT/'scripts/arb_mesh_check_forecast.py', ROOT/'scripts/mesh_problem_families.py', ROOT/'scripts/arb_cubic_reconstruction.py', ROOT/'scripts/arb_degree_selection.py', FIXTURE, REFERENCE, ROOT/'results/generated/24_depth_reference/20261029T044806546934Z/references.json', ROOT/'configs/frozen_networks_selection_fresh_2026-10-29.json', ROOT/'results/generated/35_fresh_selection_cases/20261029T111900086450Z/references.json']
    extra = ['arb_reconstruction.py','arb_shaped_remainder.py','arb_deep_network.py','arb_certification.py','arb_problem_families.py','interval_certification.py','18_continuous_reference_crosscheck.py']
    sources += [ROOT/'scripts'/name for name in extra]
    snapshot = out/'source_snapshot'; snapshot.mkdir()
    for src in sources:
        if src.suffix == '.py': (snapshot/src.name).write_bytes(src.read_bytes())
    config = dict(status='running', suite=args.suite, cases=[], targets=[str(t) for t in TARGETS], policies=list(POLICIES), repeats=3,
                  max_intervals=4096, seconds=30, source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                  python=platform.python_version(), flint=flint.__version__,
                  timing='Sequential; fresh cache per request; all attempts and selection included; imports, garbage collection, references, diagnostic probes and file writes excluded.',
                  limitations='All networks and analytic shapes have appeared in earlier work; evaluation parameters were reserved from the immediate development run but are not new independent network training. Forecast uses only completed certificate components, exact goal remains mandatory. Intermediate planned grids are not certified. Both targets are certificate enclosure ratios, not solution error.')
    def save(name, data):
        p=out/name; p.write_text(json.dumps(data,indent=2)+'\n'); assert json.loads(p.read_text()) == data
    def checkpoint(): save('config.json',config);save('runs.json',runs);save('checks.json',checks)
    fresh=ROOT/'configs/frozen_networks_selection_fresh_2026-10-29.json'
    if args.suite == 'development':
        fixture=FIXTURE; labels=('depth1_k1_seed7','depth3_k5_seed7','depth3_k5_seed17')
        terms=[dict(center='0.733',width='0.013',weight='1')]
        controls=[('narrow_peak',lambda:BumpSix('0.005')),('global_sine_50',lambda:SineSix(50)),('shift_0733',lambda:GaussianMixture(terms))]
    else:
        fixture=fresh; labels=('depth1_k1_seed37','depth3_k1_seed37','depth3_k5_seed37')
        shapes={
            'boundary_097':[dict(center='0.97',width='0.004',weight='1')],
            'shift_0137':[dict(center='0.137',width='0.007',weight='1')],
            'triple_015_047_083':[dict(center='0.15',width='0.008',weight='1/3'),dict(center='0.47',width='0.011',weight='1/3'),dict(center='0.83',width='0.006',weight='1/3')],
        }
        controls=[(label,lambda terms=terms:GaussianMixture(terms)) for label,terms in shapes.items()]
        config['shapes']=shapes
    sha=hashlib.sha256(fixture.read_bytes()).hexdigest()
    cases=[(label,lambda label=label: SixthOrderDeep(fixture,label,sha)) for label in labels]+controls
    refs=json.loads(REFERENCE.read_text())
    refs.update(json.loads((ROOT/'results/generated/24_depth_reference/20261029T044806546934Z/references.json').read_text()))
    refs.update(json.loads((ROOT/'results/generated/35_fresh_selection_cases/20261029T111900086450Z/references.json').read_text()))
    runs=[]; checks=[]; diagnostics=[]
    config['cases']=[x[0] for x in cases]; checkpoint()
    try:
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
                        else:
                            if label not in ('narrow_peak','global_sine_50'): raise AssertionError('Missing independent reference')
                            lo=hi=F(1,10000)
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
