"""Outcome checks for research scheduling, cache reuse and safe failure."""
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from flint import arb
from arb_degree_selection import POLICIES, run_request, gap_met, fraction
from arb_cubic_reconstruction import SixthOrderDeep
from arb_problem_families import SmoothProblem

ROOT = Path(__file__).resolve().parents[1]


class Zero(SmoothProblem):
    def __init__(self):
        self.setup(103)
    def residual(self, x):
        self.guard(); return arb(0)
    def r_third(self, x):
        self.guard(); return arb(0)
    def r_fourth(self, x):
        self.guard(); return arb(0)
    def error(self, x):
        raise AssertionError('Known solution must not be queried')


class BadQuadratic(Zero):
    def r_third(self, x):
        return arb('nan')


class BadBoth(BadQuadratic):
    def r_fourth(self, x):
        return arb('nan')


def main():
    out = ROOT/'results/generated/34_degree_selection_checks'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True, exist_ok=False)
    checks = []
    sources = [Path(__file__), ROOT/'scripts/arb_degree_selection.py']
    snapshot = out/'source_snapshot'; snapshot.mkdir()
    for p in sources: (snapshot/p.name).write_bytes(p.read_bytes())
    config = dict(status='running', source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    def save(name, value):
        p = out/name; p.write_text(json.dumps(value, indent=2)+'\n')
        if json.loads(p.read_text()) != value: raise OSError('Readback mismatch')
    def check(name, passed):
        checks.append(dict(name=name, passed=bool(passed))); save('checks.json', checks)
        if not passed: raise AssertionError(name)
    save('config.json', config)
    try:
        for policy in POLICIES:
            r = run_request(Zero, policy)
            check('zero_success_'+policy, r['success'] and fraction(r['final']['lower']) == fraction(r['final']['upper']) == 0)
            check('zero_stops_first_'+policy, r['attempts'] == 1)
        check('zero_lower_not_success', not gap_met(dict(lower=arb(0), upper=arb(1)), F(51, 50)))
        check('exact_gap_equality', gap_met(dict(lower=arb(1), upper=arb(3)/2), F(3, 2)))
        for policy in ('probe_cubic', 'forecast_once', 'forecast_race', 'forecast_guarded'):
            r = run_request(BadQuadratic, policy, cap=16)
            check('arithmetic_fallback_'+policy, r['success'] and r['final_degree'] == 3 and r['trace'][0]['status'] == 'arithmetic_unresolved')
            r = run_request(BadBoth, policy, cap=16)
            check('invalid_not_zero_'+policy, not r['success'] and r['final'] is None and r['stop_reason'] == 'mesh_limit')
        r = run_request(Zero, 'forecast_race', seconds=1e-30)
        check('no_time_not_zero', not r['success'] and r['final'] is None and r['attempts'] == 0 and r['stop_reason'] == 'time_limit')
        fixture = ROOT/'configs/frozen_networks_depth_2026-10-29.json'
        sha = hashlib.sha256(fixture.read_bytes()).hexdigest()
        def factory():
            p = SixthOrderDeep(fixture, 'depth1_k1_seed7', sha)
            p.error = lambda x: (_ for _ in ()).throw(AssertionError('No known error'))
            return p
        r = run_request(factory, 'forecast_race', cap=16)
        check('finite_budget_unresolved', not r['success'] and r['final'] is not None and r['stop_reason'] == 'mesh_limit')
        r = run_request(factory, 'forecast_race', max_attempts=1)
        check('attempt_budget_respected', not r['success'] and r['attempts'] == 1 and r['stop_reason'] == 'attempt_limit')
        yes = run_request(factory, 'probe_cubic', share_nodes=True)
        no = run_request(factory, 'probe_cubic', share_nodes=False)
        check('shared_cache_same_certificate', yes['final'] == no['final'])
        check('shared_cache_saves_residual_calls', yes['work']['residual_evaluations'] < no['work']['residual_evaluations'])
        check('shared_cache_derivative_unchanged', yes['work']['derivative_interval_evaluations'] == no['work']['derivative_interval_evaluations'])
        for name, r in [('shared', yes), ('separate', no)]:
            check('actual_calls_'+name, r['work']['problem_forward_calls'] == r['work']['residual_evaluations']+r['work']['derivative_interval_evaluations'])
        # Remove training and observed-error metadata from the input fixture.
        original = json.loads(fixture.read_text())
        row = next(x for x in original['networks'] if x['label'] == 'depth1_k1_seed7')
        minimal = out/'minimal_fixture.json'
        minimal.write_text(json.dumps(dict(networks=[{k:row[k] for k in ('label', 'frequency', 'layers')}]))+'\n')
        mini_sha = hashlib.sha256(minimal.read_bytes()).hexdigest()
        r = run_request(lambda: SixthOrderDeep(minimal, row['label'], mini_sha), 'probe_cubic')
        check('no_training_or_reference_metadata_needed', r['final'] == yes['final'])
        for override in [dict(policy='degree4'), dict(target=1.02), dict(initial=True), dict(cap=0),
                         dict(max_attempts=0), dict(seconds=float('nan')), dict(seconds=-1)]:
            settings = dict(policy='forecast_race'); settings.update(override)
            rejected = False
            try: run_request(Zero, **settings)
            except ValueError: rejected = True
            check('reject_'+str(override), rejected)
        config.update(status='complete', checks_passed=len(checks)); save('config.json', config)
    except Exception as error:
        config.update(status='failed', failure=repr(error)); save('config.json', config); raise
    print('COMPLETED', out, 'checks', len(checks), flush=True)


if __name__ == '__main__': main()
