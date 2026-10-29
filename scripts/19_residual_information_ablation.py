"""Compare three safe Poisson bounds using identical quadratic residual data.

Global residual maximum, spatial absolute envelope, and signed reconstruction
use the same residual samples and r''' interval bounds. No runtime superiority
is claimed: this experiment isolates information loss, not end-to-end speed.
"""
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path

from flint import arb
from arb_certification import finite, exact_fraction, export
from arb_problem_families import SmoothProblem
from arb_reconstruction import Verifier, norm_enclosure

ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / 'scripts/17_reconstruction_adaptivity_study.py'
spec = importlib.util.spec_from_file_location('study17', DRIVER)
study = importlib.util.module_from_spec(spec)
spec.loader.exec_module(study)
REFERENCE = ROOT / 'results/generated/17_reconstruction_adaptivity/20261029T034336591977Z/independent_reference_20261029T034918489790Z/references.json'


def fraction(record):
    return Fraction(int(record['numerator']), int(record['denominator']))


def envelope_bounds(cells, depth=2):
    """rho_i >= sup_cell |r|. -W''=rho, W(0)=W(1)=0, |e|<=W."""
    envelopes = []
    for c in cells:
        r0, rm, r1 = c['r0'], c['rm'], c['r1']
        coeff = [r0, 2*rm-(r0+r1)/2, r1]
        _, polynomial_upper = norm_enclosure(coeff, depth)
        rho = finite(polynomial_upper+c['delta'], 'absolute residual envelope').upper()
        envelopes.append(rho)
    maximum = max(envelopes)
    global_upper = finite(maximum/8, 'maximum residual bound').upper()
    slope = arb(0)
    for c, rho in zip(cells, envelopes):
        slope = finite(slope+rho*c['h']*(1-(c['a']+c['b'])/2), 'envelope initial slope')
    value = arb(0)
    upper = arb(0)
    for c, rho in zip(cells, envelopes):
        h = c['h']
        coeff = [value, value+slope*h/2, value+slope*h-rho*h*h/2]
        _, local_upper = norm_enclosure(coeff, depth)
        upper = max(upper, local_upper)
        value = finite(coeff[-1], 'envelope value')
        slope = finite(slope-rho*h, 'envelope slope')
    if not value.contains(0):
        raise ArithmeticError('Absolute envelope endpoint misses zero')
    return global_upper, finite(upper, 'spatial envelope norm').upper()


class ConstantResidual(SmoothProblem):
    def __init__(self, zero=False):
        self.setup(103)
        self.amplitude = arb(0) if zero else arb(1)/100
    def residual(self, x):
        self.guard()
        return self.amplitude
    def r_third(self, x):
        self.guard()
        return arb(0)


def main():
    out = ROOT/'results/generated/19_residual_information'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True, exist_ok=False)
    references = json.loads(REFERENCE.read_text())
    checks, rows = [], []

    def check(name, passed):
        checks.append(dict(name=name, passed=bool(passed)))
        (out/'checks.json').write_text(json.dumps(checks, indent=2)+'\n')
        if not passed:
            raise AssertionError(name)

    for zero in (False, True):
        p = ConstantResidual(zero)
        v = Verifier(p, 2, 2)
        v.uniform(8)
        r = v.certify()
        global_upper, absolute_upper = envelope_bounds(v.prepare(v.nodes))
        exact = Fraction(0) if zero else Fraction(1, 800)
        check(f'constant_exact_range_{zero}', exact_fraction(r['lower']) <= exact <= exact_fraction(r['upper']))
        for name, upper in [('global', global_upper), ('absolute', absolute_upper), ('signed', r['upper'])]:
            gap = exact_fraction(upper)-exact
            check(f'constant_saturation_{zero}_{name}', 0 <= gap < Fraction(1, 10**25))

    for name, factory in study.factories().items():
        for budget in (129, 513, 1025):
            p = factory()
            v = Verifier(p, 2, 2)
            v.uniform((budget-1)//2)
            r = v.certify()
            counts = (v.residual_evaluations, v.curvature_evaluations)
            global_upper, absolute_upper = envelope_bounds(v.prepare(v.nodes))
            check(f'no_extra_model_calls_{name}_{budget}', counts == (v.residual_evaluations, v.curvature_evaluations))
            check(f'exact_sample_budget_{name}_{budget}', counts[0] == budget)
            if hasattr(p, 'known_error'):
                reference_lower = reference_upper = p.known_error
                reference_kind = 'analytic_exact_maximum'
            else:
                reference_lower = fraction(references[name]['lower'])
                reference_upper = fraction(references[name]['upper'])
                reference_kind = 'independent_mpmath_continuous_enclosure'
            check(f'signed_lower_{name}_{budget}', exact_fraction(r['lower']) <= reference_lower)
            uppers = dict(global_residual_upper=global_upper, spatial_absolute_upper=absolute_upper, signed_upper=r['upper'])
            for key, upper in uppers.items():
                check(f'reference_contained_{name}_{budget}_{key}', exact_fraction(upper) >= reference_upper)
            row = dict(problem=name, sample_budget=budget, intervals=r['intervals'], degree=2, subdivision=2,
                       residual_evaluations=counts[0], derivative_interval_evaluations=counts[1],
                       signed_lower=export(r['lower'], 'down'), reference_kind=reference_kind,
                       reference_lower=str(reference_lower), reference_upper=str(reference_upper),
                       **{key: export(value, 'up') for key, value in uppers.items()})
            row['diagnostic_upper_over_reference_midpoint'] = {
                key: float(exact_fraction(value)/((reference_lower+reference_upper)/2))
                for key, value in uppers.items()}
            rows.append(row)
        (out/'rows.json').write_text(json.dumps(rows, indent=2)+'\n')
        last = rows[-1]['diagnostic_upper_over_reference_midpoint']
        print(name, '1025 samples:', {key: round(value, 6) for key, value in last.items()}, flush=True)
    sources = [Path(__file__), DRIVER, ROOT/'scripts/arb_problem_families.py', ROOT/'scripts/arb_reconstruction.py', ROOT/'scripts/arb_certification.py']
    files = sources+[ROOT/'configs/frozen_networks_round07.json', REFERENCE]
    config = dict(configurations=len(rows), checks_passed=len(checks), bits=103,
                  protocol='Same 2N+1 residual samples and N residual-third-derivative bounds; both baseline postprocessings reuse cached cells; no timing comparison.',
                  scope='Zero-Dirichlet unit-interval Poisson; smooth candidate; conditional on interval backend correctness; no independent/formal software audit.',
                  reference_use='Validation and diagnostic ratios only; never bound construction or mesh selection.',
                  sources={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    (out/'config.json').write_text(json.dumps(config, indent=2)+'\n')
    print('COMPLETED', out, 'checks', len(checks), flush=True)


if __name__ == '__main__':
    main()
