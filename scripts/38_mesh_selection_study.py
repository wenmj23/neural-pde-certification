"""Matched-target starting-mesh/local-refinement study, with all costs included."""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as F
import gc
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import statistics

import flint
from arb_mesh_selection import POLICIES, request
from arb_degree_selection import fraction
from mesh_problem_families import GaussianMixture, MPGaussianReference

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('mesh_cases33', ROOT/'scripts/33_degree_selection_study.py')
cases33 = importlib.util.module_from_spec(spec); spec.loader.exec_module(cases33)
DEFAULT_POLICIES = ['cubic_uniform16', 'cubic_uniform1', 'cubic_adaptive16', 'cubic_adaptive1',
                    'quadratic_adaptive1', 'cubic_heap16', 'cubic_heap1', 'auto_heap16', 'auto_heap1']


def shape_terms(suite):
    if suite == 'development':
        return {
            'shift_017': [dict(center='0.17', width='0.01', weight='1')],
            'double_023_079': [dict(center='0.23', width='0.009', weight='3/5'),
                               dict(center='0.79', width='0.016', weight='2/5')],
        }
    return {
        'shift_0137': [dict(center='0.137', width='0.007', weight='1')],
        'shift_0733': [dict(center='0.733', width='0.013', weight='1')],
        'boundary_097': [dict(center='0.97', width='0.004', weight='1')],
        'double_012_088': [dict(center='0.12', width='0.006', weight='1/2'),
                           dict(center='0.88', width='0.009', weight='1/2')],
        'triple_015_047_083': [dict(center='0.15', width='0.008', weight='1/3'),
                              dict(center='0.47', width='0.011', weight='1/3'),
                              dict(center='0.83', width='0.006', weight='1/3')],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite', choices=('development', 'evaluation'), default='development')
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--policies', nargs='+', choices=POLICIES, default=DEFAULT_POLICIES)
    args = parser.parse_args()
    if not 1 <= args.repeats <= 9 or len(set(args.policies)) != len(args.policies): raise ValueError('Invalid repetition/policy list')
    out = ROOT/'results/generated/38_mesh_selection'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True, exist_ok=False)
    shapes = shape_terms(args.suite)
    tasks = cases33.cases(args.suite)
    fixture, ref_path = cases33.paths(args.suite)
    references = json.loads(ref_path.read_text())
    inputs = [fixture, ref_path]
    if args.suite == 'evaluation':
        fresh_fixture, fresh_ref = cases33.paths('fresh')
        tasks += cases33.cases('fresh', fresh_fixture)
        references.update(json.loads(fresh_ref.read_text())); inputs += [fresh_fixture, fresh_ref]
    for label, terms in shapes.items():
        tasks.append((label, lambda terms=terms: GaussianMixture(terms), F(51, 50)))
    names = ['arb_mesh_selection.py', 'mesh_problem_families.py', 'arb_degree_selection.py',
             'arb_cubic_reconstruction.py', 'arb_shaped_remainder.py', 'arb_reconstruction.py',
             'arb_deep_network.py', 'arb_certification.py', 'arb_problem_families.py',
             'interval_certification.py', '18_continuous_reference_crosscheck.py', '33_degree_selection_study.py']
    sources = [Path(__file__), *[ROOT/'scripts'/name for name in names], *inputs]
    snapshot = out/'source_snapshot'; snapshot.mkdir()
    for p in sources:
        if p.suffix == '.py': (snapshot/p.name).write_bytes(p.read_bytes())
    config = dict(status='running', suite=args.suite, repeats=args.repeats, policies=args.policies,
                  cases=[dict(label=k, target=str(t)) for k, _, t in tasks], manufactured_shapes=shapes,
                  limits=dict(cap=4096, max_attempts=40, seconds=30), bits=103,
                  timing='Fresh oracle/cache, every certificate and endpoint export, diagnostics, all refinement, priority queues and startup probes included. Imports, pre-request GC, reference computation/checks and file writes excluded.',
                  selection='Cached remainder contributions only. Fixed threshold effective-cell fraction<1/3 at N>=16; auto_heap may accept locality at N=4 with threshold0.6. Small-start auto uses uniform1->4 until selection. No true error/frequency/label used.',
                  split='Development seed7+7 historical analytic cases+2 shifted shapes. Evaluation old seeds17/27, six previously frozen seed37 cases, two historical analytic parameters and five new shifted/multipeak shapes. Same PDE/network family; not new independent network training.',
                  scope='Degree at most3, smooth zero-boundary one-dimensional setting; conditional backend inclusion. Original validated kernels unchanged; heap child path checked against original on rational meshes.',
                  versions=dict(python=platform.python_version(), python_flint=flint.__version__, flint=flint.__FLINT_VERSION__),
                  source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    def save(name, value):
        p = out/name; p.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
        if json.loads(p.read_text(encoding='utf-8')) != value: raise OSError('Readback mismatch: '+name)
    save('config.json', config)
    checks, runs = [], []
    def check(name, passed):
        checks.append(dict(name=name, passed=bool(passed))); save('checks.json', checks)
        if not passed: raise AssertionError(name)
    try:
        # Known manufactured errors are confined to this reference-only path.
        for label, terms in shapes.items():
            ref = MPGaussianReference(terms).maximum(relative=F(1, 10**8), absolute=F(1, 10**14), cap=5000)
            references[label] = ref
            check('independent_shift_reference_'+label, ref['converged'])
            save('references.json', references)
        save('references.json', references)
        for repeat in range(args.repeats):
            # Rotate methods to reduce a fixed ordering advantage. Preserve all
            # observations; three repeats are not confidence intervals.
            shift = (3*repeat) % len(args.policies)
            order = args.policies[shift:]+args.policies[:shift]
            for label, factory, target in tasks:
                for policy in order:
                    gc.collect(); result = request(factory, policy, target)
                    result.update(label=label, repeat=repeat)
                    ref = references.get(label)
                    lo_ref, hi_ref = (fraction(ref['lower']), fraction(ref['upper'])) if ref else (F(1, 10000), F(1, 10000))
                    final = result['final']
                    if final:
                        lo, hi = fraction(final['lower']), fraction(final['upper'])
                        check(f'continuous_reference_{repeat}_{label}_{policy}', lo <= lo_ref and hi >= hi_ref)
                        check(f'exact_goal_{repeat}_{label}_{policy}', result['success'] == (lo == hi == 0 or (lo > 0 and hi <= target*lo)))
                        mesh = [F(x) for x in result['final_nodes']]
                        check(f'complete_mesh_{repeat}_{label}_{policy}', len(mesh) == final['intervals']+1 and mesh[0] == 0 and mesh[-1] == 1 and all(a < b for a, b in zip(mesh[:-1], mesh[1:])))
                    calls = result['work']['problem_forward_calls']
                    if calls is not None:
                        check(f'actual_calls_{repeat}_{label}_{policy}', calls == result['work']['residual_evaluations']+result['work']['derivative_interval_evaluations'])
                    runs.append(result); save('runs.json', runs)
                    print('MESH', args.suite, repeat, label, policy, result['stop_reason'],
                          'N', final['intervals'] if final else None, 'mode', result['locked_mesh'],
                          'ms', round(result['wall_ms'], 2), flush=True)
        summary = []
        for label, _, target in tasks:
            for policy in args.policies:
                group = [r for r in runs if r['label'] == label and r['policy'] == policy]
                times = [r['wall_ms'] for r in group]
                summary.append(dict(label=label, policy=policy, target=str(target),
                                    successes=sum(r['success'] for r in group), repeats=len(group),
                                    median_wall_ms=statistics.median(times), min_wall_ms=min(times), max_wall_ms=max(times),
                                    final_intervals=[r['final']['intervals'] if r['final'] else None for r in group],
                                    selected_meshes=[r['locked_mesh'] for r in group],
                                    median_residual_calls=statistics.median(r['work']['residual_evaluations'] for r in group),
                                    median_derivative_calls=statistics.median(r['work']['derivative_interval_evaluations'] for r in group)))
        save('summary.json', summary)
        for name, sha in config['source_sha256'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != sha: raise RuntimeError('Source/input changed during study: '+name)
        config.update(status='complete', requests=len(runs), successes=sum(r['success'] for r in runs), checks_passed=len(checks))
        save('config.json', config)
    except Exception as error:
        config.update(status='failed', failure=repr(error), requests_saved=len(runs), checks_saved=len(checks)); save('config.json', config); raise
    print('COMPLETED', out, 'successes', config['successes'], '/', len(runs), 'checks', len(checks), flush=True)


if __name__ == '__main__': main()
