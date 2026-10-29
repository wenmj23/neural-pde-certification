"""Measured auto-degree pilots against fixed and simple-switch controls.

No public-CLI changes. File writes, GC before requests and validation-only
references are outside the timer; all pilots, choices, certification, cache
construction and endpoint exports are inside. CPU benchmarks run sequentially.
"""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as F
import gc
import hashlib
import json
from pathlib import Path
import platform
import statistics

import flint
from arb_degree_selection import POLICIES, run_request, fraction
from arb_cubic_reconstruction import SixthOrderDeep, SineSix, BumpSix
from arb_problem_families import QuarticError
from flint import arb

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT/'configs/frozen_networks_depth_2026-10-29.json'
REFERENCE = ROOT/'results/generated/24_depth_reference/20261029T044806546934Z/references.json'


class QuarticSix(QuarticError):
    def r_fourth(self, x):
        self.guard()
        return arb(0)


def paths(suite):
    if suite != 'fresh':
        return FIXTURE, REFERENCE
    fixture = ROOT/'configs/frozen_networks_selection_fresh_2026-10-29.json'
    data = json.loads(fixture.read_text())
    for name, expected in data['config']['locked_selection_sources'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Selection code changed after fresh-case lock: '+name)
    run = ROOT/data['config']['reference_run']
    if json.loads((run/'config.json').read_text())['status'] != 'complete':
        raise ValueError('Fresh references incomplete')
    return fixture, run/'references.json'


def cases(suite, fixture=FIXTURE):
    sha = hashlib.sha256(fixture.read_bytes()).hexdigest()
    if suite == 'fresh':
        return [(row['label'], lambda label=row['label']: SixthOrderDeep(fixture, label, sha), F(51, 50))
                for row in json.loads(fixture.read_text())['networks']]
    seeds = (7,) if suite == 'development' else (17, 27)
    items = []
    for depth in (1, 2, 3):
        for frequency in (1, 5):
            for seed in seeds:
                label = f'depth{depth}_k{frequency}_seed{seed}'
                items.append((label, lambda label=label: SixthOrderDeep(fixture, label, sha), F(51, 50)))
    if suite == 'development':
        items += [('sine_1', lambda: SineSix(1), F(51, 50)),
                  ('sine_20', lambda: SineSix(20), F(6, 5)),
                  ('sine_50', lambda: SineSix(50), F(6, 5)),
                  ('mixed', lambda: SineSix(1, frequencies=[1, 5, 81],
                                          coefficients=[F(120, 128), F(7, 128), F(1, 128)]), F(51, 50)),
                  ('bump_002', lambda: BumpSix('0.02'), F(51, 50)),
                  ('bump_0005', lambda: BumpSix('0.005'), F(51, 50)),
                  ('quartic', QuarticSix, F(51, 50))]
    else:
        items += [('sine_120', lambda: SineSix(120), F(6, 5)),
                  ('bump_001', lambda: BumpSix('0.01'), F(51, 50))]
    return items


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite', choices=('development', 'evaluation', 'fresh'), default='development')
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--policies', nargs='+', choices=POLICIES, default=list(POLICIES))
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    if not 1 <= args.repeats <= 9 or len(set(args.policies)) != len(args.policies):
        raise ValueError('Invalid repeats or duplicate policies')
    out = args.out or ROOT/'results/generated/33_degree_selection'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True, exist_ok=False)
    fixture, reference = paths(args.suite)
    sources = [Path(__file__), ROOT/'scripts/arb_degree_selection.py',
               *[ROOT/'scripts'/name for name in ('arb_cubic_reconstruction.py', 'arb_deep_network.py',
                  'arb_shaped_remainder.py', 'arb_reconstruction.py', 'arb_certification.py',
                  'arb_problem_families.py', 'interval_certification.py')], fixture, reference]
    snapshot = out/'source_snapshot'; snapshot.mkdir()
    for path in sources:
        if path.suffix == '.py':
            (snapshot/path.name).write_bytes(path.read_bytes())
    tasks = cases(args.suite, fixture)
    config = dict(status='running', suite=args.suite, repeats=args.repeats,
                  policies=args.policies, cases=[dict(label=k, target=str(t)) for k, _, t in tasks],
                  limits=dict(initial=16, cap=4096, seconds=30, max_attempts=20),
                  schedule='Default interval doubling starts at 16; microprobe_cubic tries quadratic N=1 first, *_from1 fixed baselines start N=1. Same exact gap target. Pilot/starting budgets are part of the policy.',
                  selection='Only degree 2/3; uniform meshes. Certificate Z/D data and measured per-attempt cost. No exact solution/reference used by selector.',
                  timing='Fresh oracle/cache and all pilots, decisions, failed attempts, certificate exports included. Imports, pre-request GC, file writes and external references excluded.',
                  split=('Fresh seed37 final networks generated after selector/driver sources locked; no best-output selection or post-evaluation selector retuning.' if args.suite == 'fresh' else
                         'Prespecified seed7 development, seed17/27 evaluation within previously studied fixtures; not an independent unseen network holdout. Evaluation adds two analytic parameters.'),
                  versions=dict(python=platform.python_version(), python_flint=flint.__version__, flint=flint.__FLINT_VERSION__),
                  source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    def save(name, value):
        p = out/name
        p.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
        if json.loads(p.read_text(encoding='utf-8')) != value:
            raise OSError('Readback mismatch: '+name)
    save('config.json', config)
    refs = json.loads(reference.read_text())
    runs, checks = [], []
    def check(name, passed):
        checks.append(dict(name=name, passed=bool(passed)))
        save('checks.json', checks)
        if not passed:
            raise AssertionError(name)
    try:
        for repeat in range(args.repeats):
            order = args.policies[repeat:]+args.policies[:repeat]
            for label, factory, target in tasks:
                for policy in order:
                    gc.collect()
                    result = run_request(factory, policy, target)
                    result.update(label=label, repeat=repeat)
                    final = result['final']
                    if final is not None:
                        lo, hi = fraction(final['lower']), fraction(final['upper'])
                        if label in refs:
                            lower, upper = fraction(refs[label]['lower']), fraction(refs[label]['upper'])
                        else:
                            lower = upper = F(1, 10000)
                        check(f'continuous_reference_{repeat}_{label}_{policy}', lo <= lower and hi >= upper)
                        check(f'exact_goal_{repeat}_{label}_{policy}', result['success'] == (lo == hi == 0 or (lo > 0 and hi <= target*lo)))
                    check(f'degree_cap_{repeat}_{label}_{policy}', all(row['degree'] in (2, 3) for row in result['trace']))
                    calls = result['work']['problem_forward_calls']
                    if calls is not None:
                        check(f'actual_oracle_calls_{repeat}_{label}_{policy}', calls == result['work']['residual_evaluations']+result['work']['derivative_interval_evaluations'])
                    runs.append(result); save('runs.json', runs)
                    print('SELECTOR', args.suite, repeat, label, policy, result['stop_reason'],
                          'degree', result['final_degree'], 'N', final['intervals'] if final else None,
                          'ms', round(result['wall_ms'], 2), flush=True)
        summary = []
        for label, _, target in tasks:
            medians = {}
            for policy in args.policies:
                group = [r for r in runs if r['label'] == label and r['policy'] == policy]
                times = [r['wall_ms'] for r in group]
                medians[policy] = statistics.median(times)
                summary.append(dict(label=label, target=str(target), policy=policy,
                                    successes=sum(r['success'] for r in group), repeats=len(group),
                                    median_wall_ms=statistics.median(times), min_wall_ms=min(times), max_wall_ms=max(times),
                                    final_degrees=[r['final_degree'] for r in group],
                                    final_intervals=[r['final']['intervals'] if r['final'] else None for r in group]))
            fixed = [value for policy, value in medians.items() if policy.startswith('fixed')]
            if fixed:
                best = min(fixed)
                for item in summary:
                    if item['label'] == label:
                        item['cost_over_best_fixed'] = item['median_wall_ms']/best
        save('summary.json', summary)
        config.update(status='complete', checks_passed=len(checks), timed_requests=len(runs),
                      successes=sum(r['success'] for r in runs))
        save('config.json', config)
    except Exception as error:
        config.update(status='failed', failure=repr(error), checks_saved=len(checks), requests_saved=len(runs))
        save('config.json', config)
        raise
    print('COMPLETED', out, 'successes', config['successes'], '/', len(runs), 'checks', len(checks), flush=True)


if __name__ == '__main__':
    main()
