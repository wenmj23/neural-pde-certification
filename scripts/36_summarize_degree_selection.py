"""Audit all saved scheduler stages and summarize matched comparisons.

Reference/error information is used only here, after all timed decisions.
Historical source hashes are checked against each run's source snapshots.
"""
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parents[1]


def frac(item):
    return F(int(item['numerator']), int(item['denominator']))


def main():
    runs = []
    for path in sorted((ROOT/'results/generated/33_degree_selection').glob('*/config.json')):
        config = json.loads(path.read_text())
        if config['status'] == 'complete': runs.append((path.parent, config))
    if not any(c['suite'] == 'fresh' for _, c in runs):
        raise ValueError('Fresh locked-rule validation must complete first')
    out = ROOT/'results/generated/36_degree_selection_summary'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True, exist_ok=False)
    (out/'executed_source.py').write_bytes(Path(__file__).read_bytes())
    stages, intermediate, requests, total_checks = [], 0, 0, 0
    labels = set()
    for directory, config in runs:
        reference_path = next(ROOT/name for name in config['source_sha256'] if name.endswith('/references.json'))
        refs = json.loads(reference_path.read_text())
        for name, sha in config['source_sha256'].items():
            original = ROOT/name
            p = directory/'source_snapshot'/original.name if original.suffix == '.py' else original
            if hashlib.sha256(p.read_bytes()).hexdigest() != sha:
                raise AssertionError('Source/input mismatch: '+str(p))
        results = json.loads((directory/'runs.json').read_text())
        checks = json.loads((directory/'checks.json').read_text())
        assert len(results) == config['timed_requests']
        assert len(checks) == config['checks_passed'] and all(x['passed'] for x in checks)
        assert sum(r['success'] for r in results) == config['successes']
        for r in results:
            labels.add(r['label'])
            ref = refs.get(r['label'])
            reference_lo, reference_hi = (frac(ref['lower']), frac(ref['upper'])) if ref else (F(1, 10000), F(1, 10000))
            for row in r['trace']:
                assert row['degree'] in (2, 3)
                if row['status'] != 'finite': continue
                enclosure = row['enclosure']; lo, hi = frac(enclosure['lower']), frac(enclosure['upper'])
                assert lo <= reference_lo and hi >= reference_hi
                actual = lo == hi == 0 or (lo > 0 and hi <= F(r['target'])*lo)
                assert actual == row['goal_met']
                intermediate += 1
        summary = json.loads((directory/'summary.json').read_text())
        costs = {}
        for policy in config['policies']:
            group = [x for x in summary if x['policy'] == policy]
            ratios = [x['cost_over_best_fixed'] for x in group]
            costs[policy] = dict(median_cost_over_best_fixed=statistics.median(ratios),
                                 maximum_cost_over_best_fixed=max(ratios),
                                 selected_degrees=sorted({d for x in group for d in x['final_degrees']}),
                                 median_total_ms=sum(x['median_wall_ms'] for x in group))
        stages.append(dict(directory=str(directory.relative_to(ROOT)), suite=config['suite'],
                           case_count=len(config['cases']), requests=len(results), checks=len(checks),
                           successes=config['successes'], policies=costs,
                           source_sha256=config['source_sha256']))
        requests += len(results); total_checks += len(checks)
    fresh_fixture = json.loads((ROOT/'configs/frozen_networks_selection_fresh_2026-10-29.json').read_text())
    for name, sha in fresh_fixture['config']['locked_selection_sources'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == sha
    check_path = max((ROOT/'results/generated/34_degree_selection_checks').glob('*/config.json'))
    check_config = json.loads(check_path.read_text())
    assert check_config['status'] == 'complete'
    for name, sha in check_config['source_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == sha
    pipeline = json.loads((ROOT/'results/generated/32_pipeline_checks/20261029T070210333477Z/config.json').read_text())
    for name, sha in pipeline['sources'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == sha
    digest = dict(status='complete', stages=stages, unique_case_count=len(labels),
                  timed_requests=requests, successful_requests=sum(x['successes'] for x in stages),
                  recorded_comparison_checks=total_checks,
                  intermediate_continuous_enclosures_checked=intermediate,
                  latest_scheduler_checks=check_config['checks_passed'],
                  public_pipeline_sources_unchanged=True, fresh_selection_source_lock_verified=True,
                  interpretation='Different stage protocols and repeated cases; not a single full-factorial dataset. Ratios compare only within the same stage/case. Tiny variations are not confidence intervals.',
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    p = out/'digest.json'; p.write_text(json.dumps(digest, indent=2)+'\n')
    assert json.loads(p.read_text()) == digest
    print('COMPLETED', out, 'cases', len(labels), 'requests', requests, 'checks', total_checks,
          'intermediate enclosures', intermediate, flush=True)
    for stage in stages:
        print(stage['suite'], stage['directory'])
        for policy, item in stage['policies'].items():
            print(' ', policy, 'median relative cost', round(item['median_cost_over_best_fixed'], 3),
                  'max', round(item['maximum_cost_over_best_fixed'], 3))


if __name__ == '__main__': main()
