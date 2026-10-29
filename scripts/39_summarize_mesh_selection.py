"""Audit stored meshes, continuous references and original-kernel equivalence."""
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import statistics

from arb_degree_selection import pack, fraction
from arb_cubic_reconstruction import CubicVerifier
from arb_shaped_remainder import ShapedVerifier

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('mesh_study38', ROOT/'scripts/38_mesh_selection_study.py')
study = importlib.util.module_from_spec(spec); spec.loader.exec_module(study)


def factories(suite):
    items = study.cases33.cases(suite)
    if suite == 'evaluation':
        p, _ = study.cases33.paths('fresh'); items += study.cases33.cases('fresh', p)
    for label, terms in study.shape_terms(suite).items():
        items.append((label, lambda terms=terms: study.GaussianMixture(terms), F(51, 50)))
    return {label:factory for label, factory, _ in items}


def main():
    studies = []
    for path in sorted((ROOT/'results/generated/38_mesh_selection').glob('*/config.json')):
        config = json.loads(path.read_text())
        if config['status'] == 'complete': studies.append((path.parent, config))
    if {c['suite'] for _, c in studies} != {'development', 'evaluation'}:
        raise ValueError('Both formal suites must complete')
    out = ROOT/'results/generated/39_mesh_summary'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True, exist_ok=False)
    (out/'executed_source.py').write_bytes(Path(__file__).read_bytes())
    stages = []; labels = set(); requests = successful = checks = intermediate = 0
    matches = []; verified = {}; selected_rows = []
    def save(name, value):
        p = out/name; p.write_text(json.dumps(value, indent=2)+'\n')
        if json.loads(p.read_text()) != value: raise OSError('Readback mismatch')
    for directory, config in studies:
        for name, expected in config['source_sha256'].items():
            original = ROOT/name
            p = directory/'source_snapshot'/original.name if original.suffix == '.py' else original
            assert hashlib.sha256(p.read_bytes()).hexdigest() == expected, name
        records = json.loads((directory/'runs.json').read_text())
        saved_checks = json.loads((directory/'checks.json').read_text())
        assert len(records) == config['requests']
        assert len(saved_checks) == config['checks_passed'] and all(x['passed'] for x in saved_checks)
        assert sum(r['success'] for r in records) == config['successes']
        refs = json.loads((directory/'references.json').read_text())
        for r in records:
            labels.add(r['label']); ref = refs.get(r['label'])
            reflo, refhi = (fraction(ref['lower']), fraction(ref['upper'])) if ref else (F(1, 10000), F(1, 10000))
            for row in r['trace']:
                if row['status'] != 'finite': continue
                lo, hi = fraction(row['enclosure']['lower']), fraction(row['enclosure']['upper'])
                assert lo <= reflo and hi >= refhi, (r['label'], r['policy'], row['attempt'])
                actual = lo == hi == 0 or (lo > 0 and hi <= F(r['target'])*lo)
                assert row['goal_met'] == actual
                intermediate += 1
        summary = json.loads((directory/'summary.json').read_text())
        grouped = {}
        for row in summary: grouped.setdefault(row['label'], {})[row['policy']] = row
        ratios = {}
        for label, group in grouped.items():
            fixed = [x['median_wall_ms'] for name, x in group.items()
                     if not name.startswith('auto') and x['successes'] == x['repeats']]
            if not fixed: raise AssertionError('No successful fixed baseline for '+label)
            best = min(fixed)
            for name, row in group.items():
                ratio = row['median_wall_ms']/best if row['successes'] == row['repeats'] else None
                if ratio is not None: ratios.setdefault(name, []).append(ratio)
                selected_rows.append(dict(suite=config['suite'], label=label, policy=name,
                                          median_wall_ms=row['median_wall_ms'], min_wall_ms=row['min_wall_ms'], max_wall_ms=row['max_wall_ms'],
                                          relative_to_best_fixed=ratio, final_intervals=row['final_intervals'],
                                          selected_meshes=row['selected_meshes']))
        stats = {name:dict(median_cost_over_best_fixed=statistics.median(values), max_cost_over_best_fixed=max(values))
                 for name, values in ratios.items()}
        stages.append(dict(directory=str(directory.relative_to(ROOT)), suite=config['suite'], cases=len(config['cases']),
                           requests=len(records), checks=len(saved_checks), successes=config['successes'], relative_costs=stats))
        requests += len(records); checks += len(saved_checks); successful += config['successes']
        # All heap-produced final grids are recomputed by the original kernels.
        # Repeated identical meshes are checked once, with final endpoints of
        # every repeated result checked against that same recomputed enclosure.
        fs = factories(config['suite'])
        for r in records:
            if 'heap' not in r['policy'] or r['final'] is None: continue
            mesh_hash = hashlib.sha256(json.dumps(r['final_nodes']).encode()).hexdigest()
            key = (r['label'], r['degree'], mesh_hash)
            if key not in verified:
                p = fs[r['label']](); v = CubicVerifier(p, 'moment') if r['degree'] == 3 else ShapedVerifier(p, 'moment')
                recomputed = pack(v.certify([F(x) for x in r['final_nodes']]))
                verified[key] = recomputed
                matches.append(dict(label=r['label'], degree=r['degree'], mesh_sha256=mesh_hash, matched=True))
                print('ORIGINAL_MESH_CHECKED', r['label'], r['policy'], r['final']['intervals'], flush=True)
            assert verified[key] == r['final'], (r['label'], r['policy'])
            save('original_kernel_matches.json', matches)
    latest = max((ROOT/'results/generated/37_mesh_selection_checks').glob('*/config.json'))
    procedure = json.loads(latest.read_text()); assert procedure['status'] == 'complete'
    for name, sha in procedure['source_sha256'].items(): assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == sha
    old_pipeline = json.loads((ROOT/'results/generated/32_pipeline_checks/20261029T070210333477Z/config.json').read_text())
    for name, sha in old_pipeline['sources'].items(): assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == sha
    result = dict(status='complete', stages=stages, unique_cases=len(labels), requests=requests, successful_requests=successful,
                  comparison_checks=checks, intermediate_ranges_checked=intermediate, original_final_meshes_checked=len(matches),
                  latest_procedure_checks=procedure['checks_passed'], public_cli_sources_unchanged=True,
                  note='Within-stage matched-target ratios; three-repeat medians/min/max, no confidence intervals. Models include older fixtures; synthetic peaks are not PINNs. All guarantees conditional on backend inclusion.',
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    save('digest.json', result); save('comparisons.json', selected_rows)
    print('COMPLETED', out, 'requests', requests, 'success', successful, 'cases', len(labels), 'checks', checks,
          'intermediate', intermediate, 'original meshes', len(matches), flush=True)


if __name__ == '__main__': main()
