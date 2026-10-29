"""Check new heap cells against original kernels and manufactured references."""
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import mpmath as mp
from flint import arb

from arb_mesh_selection import POLICIES, RecordedCubic, indicators, bulk_nodes, request
from arb_cubic_reconstruction import CubicVerifier, SineSix, BumpSix, SixthOrderDeep
from arb_degree_selection import pack, fraction
from mesh_problem_families import GaussianMixture, MPGaussianReference

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('checks34', ROOT/'scripts/34_degree_selection_checks.py')
old = importlib.util.module_from_spec(spec); spec.loader.exec_module(old)


def main():
    out = ROOT/'results/generated/37_mesh_selection_checks'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True, exist_ok=False)
    sources = [Path(__file__), *[ROOT/'scripts'/n for n in ('arb_mesh_selection.py', 'mesh_problem_families.py',
               'arb_cubic_reconstruction.py', 'arb_reconstruction.py', 'arb_shaped_remainder.py',
               'arb_deep_network.py', '34_degree_selection_checks.py')]]
    snapshot = out/'source_snapshot'; snapshot.mkdir()
    for p in sources: (snapshot/p.name).write_bytes(p.read_bytes())
    config = dict(status='running', source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    checks = []
    def save(name, value):
        p = out/name; p.write_text(json.dumps(value, indent=2)+'\n')
        if json.loads(p.read_text()) != value: raise OSError('Readback mismatch')
    def check(name, value):
        checks.append(dict(name=name, passed=bool(value))); save('checks.json', checks)
        if not value: raise AssertionError(name)
    save('config.json', config)
    try:
        fixture = ROOT/'configs/frozen_networks_depth_2026-10-29.json'
        sha = hashlib.sha256(fixture.read_bytes()).hexdigest()
        factories = [('sine', lambda: SineSix(7)), ('narrow', lambda: BumpSix('0.005'))]
        for depth in (1, 2, 3):
            label = f'depth{depth}_k1_seed7'
            factories.append((label, lambda label=label: SixthOrderDeep(fixture, label, sha)))
        meshes = [[F(0), F(1, 8), F(3, 16), F(1, 2), F(7, 8), F(1)],
                  [F(0), F(1, 7), F(2, 5), F(4, 5), F(1)]]
        for label, factory in factories:
            for index, mesh in enumerate(meshes):
                oracle = factory(); v = RecordedCubic(oracle)
                for a, b in zip(mesh[:-1], mesh[1:]): v.cell(a, b)
                counts = (v.residual_evaluations, v.curvature_evaluations)
                v.nodes = mesh; result = v.certify()
                check(f'cached_cell_no_extra_calls_{label}_{index}', counts == (v.residual_evaluations, v.curvature_evaluations))
                reference = CubicVerifier(factory()).certify(mesh)
                check(f'child_builder_matches_original_{label}_{index}', pack(result) == pack(reference))
                counts = (v.residual_evaluations, v.curvature_evaluations)
                scores, concentration = indicators(v)
                check(f'finite_scaled_scores_{label}_{index}', all(math.isfinite(x) and 0 <= x <= 1 for x in scores) and 0 < concentration <= 1)
                check(f'indicators_no_oracle_calls_{label}_{index}', counts == (v.residual_evaluations, v.curvature_evaluations))
            v = RecordedCubic(factory()); v.uniform(4); v.certify(); v.heap_to(32)
            check('heap_valid_mesh_'+label, len(v.nodes) == 33 and v.nodes[0] == 0 and v.nodes[-1] == 1 and all(a < b for a, b in zip(v.nodes[:-1], v.nodes[1:])))
            rejected = False
            try: indicators(v)
            except ValueError: rejected = True
            check('stale_indicator_rejected_'+label, rejected)
            result = v.certify(); reference = CubicVerifier(factory()).certify(v.nodes)
            check('heap_mesh_original_certificate_'+label, pack(result) == pack(reference))
        for policy in POLICIES:
            r = request(old.Zero, policy)
            check('zero_'+policy, r['success'] and fraction(r['final']['lower']) == fraction(r['final']['upper']) == 0 and r['attempts'] == 1)
            r = request(old.BadBoth, policy, cap=16)
            check('invalid_not_zero_'+policy, not r['success'] and r['final'] is None and r['stop_reason'] == 'mesh_limit')
        r = request(old.Zero, 'auto_heap1', seconds=1e-30)
        check('no_time_no_fake_zero', not r['success'] and r['final'] is None and r['attempts'] == 0)
        r = request(lambda: BumpSix('0.005'), 'auto_heap16', max_attempts=1)
        check('attempt_limit_keeps_completed_mesh', not r['success'] and r['attempts'] == 1 and len(r['final_nodes']) == 17 and r['final']['intervals'] == 16)
        r = request(lambda: BumpSix('0.005'), 'auto_heap16', cap=16)
        check('mesh_limit_not_rejection', not r['success'] and r['final'] is not None and r['stop_reason'] == 'mesh_limit')
        # Independent off-center/multi-peak checks; true error used only here.
        terms = [dict(center='0.17', width='0.01', weight='3/5'),
                 dict(center='0.81', width='0.02', weight='2/5')]
        mp.mp.dps = 140
        problem = GaussianMixture(terms); function = problem.mp_function()
        for term in terms:
            x = F(term['center'])+F(term['width'])/5
            xx = mp.mpf(x.numerator)/x.denominator; aa = arb(x.numerator)/x.denominator
            for name, order in [('residual', 2), ('r_third', 5), ('r_fourth', 6)]:
                value = -mp.diff(function, xx, order)
                check(f'shifted_derivative_{term["center"]}_{name}', getattr(problem, name)(aa).overlaps(arb(str(value))))
            check('shifted_error_value_'+term['center'], problem.error(aa).overlaps(arb(str(function(xx)))))
        reference = MPGaussianReference(terms).maximum(relative=F(1, 10**8), absolute=F(1, 10**14), cap=5000)
        check('independent_shifted_reference_converged', reference['converged'])
        for policy in ('cubic_uniform16', 'cubic_adaptive1', 'cubic_heap1', 'auto_heap1'):
            r = request(lambda: GaussianMixture(terms), policy)
            check('shifted_continuous_range_'+policy, r['success'] and fraction(r['final']['lower']) <= fraction(reference['lower']) and fraction(r['final']['upper']) >= fraction(reference['upper']))
        save('shifted_reference.json', reference)
        for score in ([-1.], [float('nan')], [float('inf')]):
            rejected = False
            try: bulk_nodes([F(0), F(1)], score, 16)
            except ValueError: rejected = True
            check('invalid_marking_'+str(score), rejected)
        config.update(status='complete', checks_passed=len(checks)); save('config.json', config)
    except Exception as error:
        config.update(status='failed', failure=repr(error), saved_checks=len(checks)); save('config.json', config); raise
    print('COMPLETED', out, 'checks', len(checks), flush=True)


if __name__ == '__main__': main()
