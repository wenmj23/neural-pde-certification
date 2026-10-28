"""Audit the scoped interval implementation with the canonical frozen PINN.

This is an assistant code/math review plus executable checks, not independent
software verification or a proof of mpmath. All runs preserve earlier evidence.
"""
import argparse
from datetime import datetime, timezone
from fractions import Fraction
import importlib.util
import inspect
import json
from pathlib import Path
import platform
import shutil
import sys
import tempfile

import mpmath as mp
from mpmath.libmp import to_rational
from interval_certification import (context, endpoint_fraction, export_endpoint,
    export_bounds, FrozenNeuralPoisson, certify, continuous_reference, sha256)

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / 'cases/poisson_1d_sine_k1_seed7'


class AnalyticProblem:
    boundary_error = (0, 0)
    def __init__(self, ctx, name):
        self.ctx, self.name = ctx, name
        self.A = ctx.mpf(1) / 10000
    def residual(self, x):
        c = self.ctx
        if self.name == 'zero': return c.mpf(0)
        if self.name == 'constant': return c.mpf(2)
        if self.name == 'linear': return 6*x
        if self.name == 'quartic': return self.A*(-32+192*x-192*x**2)
        k = int(self.name.split('_')[1])
        return self.A*(k*c.pi)**2*c.sin(k*c.pi*x)
    def r_second(self, x):
        c = self.ctx
        if self.name in ('zero', 'constant', 'linear'): return c.mpf(0)
        if self.name == 'quartic': return -384*self.A
        k = int(self.name.split('_')[1])
        return -self.A*(k*c.pi)**4*c.sin(k*c.pi*x)
    def exact_norm_enclosure(self):
        if self.name == 'zero': return self.ctx.mpf(0)
        if self.name == 'constant': return self.ctx.mpf(1)/4
        if self.name == 'linear': return 2/(3*self.ctx.sqrt(3))
        return self.A


def serial_result(r):
    out = export_bounds(r)
    for k in ('intervals', 'precision_digits', 'elapsed_s', 'method', 'scope'):
        if k in r: out[k] = r[k]
    if 'z_lower' in r:
        out['auxiliary_norm_lower'] = export_endpoint(r['z_lower'], 'down')
        out['auxiliary_norm_upper'] = export_endpoint(r['z_upper'], 'up')
        out['remainder_norm_upper'] = export_endpoint(r['remainder_upper'], 'up')
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    out = args.out or ROOT / 'results/generated/12_interval_audit' / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True, exist_ok=False)
    checks = []
    def check(name, condition, **details):
        checks.append(dict(name=name, passed=bool(condition), **details))
        (out/'checks.json').write_text(json.dumps(checks, indent=2)+'\n')
        if not condition: raise AssertionError(name)
    ctx = context(30)
    problem = FrozenNeuralPoisson(CASE, ctx)
    check('49_binary64_weights_imported_exactly', problem.input_conversions == 49)
    check('hard_zero_boundary_values', problem.solution(ctx.mpf(0)) == ctx.mpf(0)
          and problem.solution(ctx.mpf(1)) == ctx.mpf(0))
    for dps in (10, 201, 30.5, True):
        try: context(dps)
        except ValueError: rejected = True
        else: rejected = False
        check('reject_precision_'+str(dps), rejected)
    for count in (0, -1, 2.5, True):
        try: certify(problem, count)
        except ValueError: rejected = True
        else: rejected = False
        check('reject_mesh_'+str(count), rejected)
    with tempfile.TemporaryDirectory() as temp:
        directory = Path(temp)
        for name in ('source', 'boundary', 'architecture', 'tampered_weights', 'nonfinite_weight'):
            shutil.copytree(CASE, directory/name)
            target = directory/name
            case = json.loads((target/'case.json').read_text())
            if name == 'source': case['problem']['source'] = '2'
            elif name == 'boundary': case['problem']['boundary']['values'][0] = 1
            elif name == 'architecture': case['network']['architecture'] = [1, 32, 1]
            else:
                weights = json.loads((target/'weights.json').read_text())
                weights['output.bias'][0] = float('nan') if name == 'nonfinite_weight' else weights['output.bias'][0]+.1
                (target/'weights.json').write_text(json.dumps(weights))
                if name == 'nonfinite_weight': case['network']['weights_sha256'] = sha256(target/'weights.json')
            (target/'case.json').write_text(json.dumps(case))
            try: FrozenNeuralPoisson(target, ctx)
            except ValueError: rejected = True
            else: rejected = False
            check('reject_case_'+name, rejected)
    # Exact arithmetic checks for endpoint export, including signs and scales.
    for name, point in [('third_low', (ctx.mpf(1)/3).a), ('third_high', (ctx.mpf(1)/3).b),
                        ('negative', -(ctx.mpf(1)/7).a), ('zero', ctx.mpf(0)),
                        ('small', ctx.mpf('1e-100').a), ('large', ctx.mpf('1e100').b)]:
        for direction in ('down', 'up'):
            e = export_endpoint(point, direction)
            rational = Fraction(int(e['numerator']), int(e['denominator']))
            decimal = Fraction(e['decimal_outward'])
            check('export_'+name+'_'+direction, rational == endpoint_fraction(point)
                  and (decimal <= rational if direction == 'down' else decimal >= rational))
    point = (ctx.mpf(1)/3).b
    old_text = str(point)[1:-1].split(',')[0].strip()
    check('old_point_text_is_not_always_an_outward_upper_bound', Fraction(old_text) < endpoint_fraction(point),
          old_text=old_text, implication='Point formatting is display, not a safe bound serialization')
    print('Input guards and exact/outward endpoint export checks passed.', flush=True)

    known_runs = []
    for name, meshes in [('zero', [1, 32]), ('constant', [1, 8, 64]),
                         ('linear', [8, 64]), ('quartic', [16, 128]),
                         ('sine_1', [64]), ('sine_20', [256, 1000]), ('sine_200', [100])]:
        p = AnalyticProblem(ctx, name)
        target = p.exact_norm_enclosure()
        for n in meshes:
            result = certify(p, n)
            check('known_norm_'+name+'_'+str(n), endpoint_fraction(result['lower']) <= endpoint_fraction(target)
                  and endpoint_fraction(result['upper']) >= endpoint_fraction(target, 1))
            known_runs.append(dict(problem=name, **serial_result(result)))
    (out/'analytic_checks.json').write_text(json.dumps(known_runs, indent=2)+'\n')
    print('Known-norm checks passed, including an intentionally aliased sine.', flush=True)

    weights = json.loads((CASE/'weights.json').read_text())
    with mp.workdps(100):
        w = [mp.mpf(float(v[0])) for v in weights['hidden.weight']]
        b = [mp.mpf(float(v)) for v in weights['hidden.bias']]
        a = [mp.mpf(float(v)) for v in weights['output.weight'][0]]
        c = mp.mpf(float(weights['output.bias'][0]))
        def v(x): return x*(1-x)*(c+sum(aa*mp.tanh(ww*x+bb) for ww, bb, aa in zip(w,b,a)))
        def contains(iv_value, value):
            approximate = Fraction(*to_rational(value._mpf_))
            return endpoint_fraction(iv_value) <= approximate <= endpoint_fraction(iv_value, 1)
        sample_rationals = [(0,1), (1,1000), (1,10), (17,64), (1,2), (37,64), (9,10), (999,1000), (1,1)]
        for numerator, denominator in sample_rationals:
            x = mp.mpf(numerator)/denominator
            xi = ctx.mpf(numerator)/denominator
            values = [('solution', problem.solution(xi), v(x)),
                      ('residual', problem.residual(xi), -mp.diff(v,x,2)-mp.pi**2*mp.sin(mp.pi*x)),
                      ('r_second', problem.r_second(xi), -mp.diff(v,x,4)+mp.pi**4*mp.sin(mp.pi*x))]
            for label, enclosure, value in values:
                check(f'high_precision_{label}_{numerator}_{denominator}', contains(enclosure, value))
        for i in (0, 3, 7):
            cell = ctx.mpf([ctx.mpf(i)/8, ctx.mpf(i+1)/8])
            bound = problem.r_second(cell)
            for j in (0,1,2):
                x = mp.mpf(2*i+j)/16
                derivative = -mp.diff(v,x,4)+mp.pi**4*mp.sin(mp.pi*x)
                check(f'whole_cell_derivative_sample_{i}_{j}', contains(bound, derivative))
    print('Independent 100-digit point/derivative diagnostics passed (not a backend proof).', flush=True)

    records, internal = [], []
    baseline = json.loads((CASE/'baseline.json').read_text())
    for dps in (30, 50):
        ctx_run = context(dps)
        p = FrozenNeuralPoisson(CASE, ctx_run)
        for n in (128, 512, 1000):
            result = certify(p, n)
            ordinary = next(r for r in baseline['certificates'] if r['method']=='linear' and r['intervals']==n)
            check(f'float_baseline_inside_interval_{dps}_{n}',
                  endpoint_fraction(result['lower']) <= Fraction.from_float(ordinary['lower'])
                  and endpoint_fraction(result['upper']) >= Fraction.from_float(ordinary['upper']))
            records.append(serial_result(result))
            internal.append(result)
            (out/'certificates.json').write_text(json.dumps(records, indent=2)+'\n')
            print('Fixed PINN',dps,'digits,',n,'intervals:',records[-1]['lower']['decimal_outward'],
                  records[-1]['upper']['decimal_outward'],f"({result['elapsed_s']:.3f}s)", flush=True)
    # Same formulas, mesh, weights and arithmetic order as old prototype 09.
    spec = importlib.util.spec_from_file_location('old09_audit', ROOT/'scripts/09_interval_prototype.py')
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    old.iv.dps = 30
    _, lo, hi = old.run(old.NeuralResidual(weights),128)
    check('refactored_matches_original_128_dps30_exact_endpoints',
          endpoint_fraction(lo) == endpoint_fraction(internal[0]['lower'])
          and endpoint_fraction(hi) == endpoint_fraction(internal[0]['upper']))

    reference_problem = FrozenNeuralPoisson(CASE, context(50))
    reference = continuous_reference(reference_problem, 2048)
    (out/'continuous_reference.json').write_text(json.dumps(serial_result(reference), indent=2)+'\n')
    for r in internal:
        check(f'continuous_reference_inside_certificate_{r["precision_digits"]}_{r["intervals"]}',
              endpoint_fraction(r['lower']) <= endpoint_fraction(reference['lower'])
              and endpoint_fraction(r['upper']) >= endpoint_fraction(reference['upper']))
    for n in (128,512,1000):
        low = next(r for r in internal if r['precision_digits']==30 and r['intervals']==n)
        high = next(r for r in internal if r['precision_digits']==50 and r['intervals']==n)
        check('precision_runs_overlap_'+str(n),
              max(endpoint_fraction(low['lower']),endpoint_fraction(high['lower']))
              <= min(endpoint_fraction(low['upper']),endpoint_fraction(high['upper'])))
    print('Continuous reference:',export_bounds(reference),f"({reference['elapsed_s']:.3f}s)", flush=True)
    result = next(r for r in internal if r['precision_digits']==50 and r['intervals']==512)
    decisions = []
    for tolerance in (Fraction(11,100000), Fraction(1,10000)):
        lower, upper = endpoint_fraction(result['lower']), endpoint_fraction(result['upper'])
        status = 'accept' if upper <= tolerance else ('reject' if lower > tolerance else 'uncertain')
        decisions.append(dict(tolerance_numerator=str(tolerance.numerator),
                              tolerance_denominator=str(tolerance.denominator),status=status,
                              scope='conditional interval certificate of continuous real-function error'))
    check('tolerance_decisions_accept_then_reject', [d['status'] for d in decisions]==['accept','reject'])
    (out/'decisions.json').write_text(json.dumps(decisions, indent=2)+'\n')
    sources = [Path(__file__),ROOT/'scripts/interval_certification.py',ROOT/'scripts/09_interval_prototype.py']
    backend_files = [Path(inspect.getfile(mp.iv.__class__)), Path(inspect.getfile(mp.libmp.libmpi)),
                     Path(inspect.getfile(mp.libmp.libelefun)), Path(inspect.getfile(mp.libmp.libmpf))]
    metadata = dict(case_id=problem.case['case_id'],weights_sha256=problem.weights_sha256,
                    python=sys.version,platform=platform.platform(),mpmath=mp.__version__,
                    backend_files={str(p):sha256(p) for p in backend_files},
                    sources={str(p.relative_to(ROOT)):sha256(p) for p in sources},
                    case_sha256=sha256(CASE/'case.json'),checks_passed=len(checks),
                    interval_configurations=len(records),analytic_configurations=len(known_runs),
                    guarantee='Conditional on reviewed mathematical assumptions and correctness of the pinned interval backend.',
                    limitations=['mpmath iv is documented experimental',
                                 'Assistant review, no independent audit or formal verification',
                                 'Certificate targets the mathematical network with exact binary64 weights, not arbitrary floating inference programs',
                                 'Continuous reference shares the same interval backend; high-precision diagnostics are not proofs'],
                    changes_from_09=['Canonical frozen case and input/hash guards',
                                     'Local interval precision contexts',
                                     'Exact rational endpoint export and directed decimal display',
                                     'Separate continuous-reference enclosure and audit checks'])
    (out/'audit_summary.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print('Audit checks passed:',len(checks),'interval configurations:',len(records),flush=True)
    print('RESULT_DIRECTORY',out,flush=True)


if __name__=='__main__': main()
