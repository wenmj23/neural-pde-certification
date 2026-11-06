"""Check residual and implemented derivatives with controlled source mistakes."""
from fractions import Fraction as F
from math import factorial
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import tempfile
from flint import arb
from mpmath.libmp import to_rational
from arb_certification import exact_fraction
from arb_reconstruction import point

ROOT = Path(__file__).resolve().parents[1]
CFG = ROOT / 'configs/multiorder_fault_diagnostic_2026-11-06.json'
BASE = ROOT / 'results/2026-11-06'
OUT = BASE / 'multiorder_fault_diagnostic'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def load(path):
    spec = importlib.util.spec_from_file_location('mp_oracle_107', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def source(coeff, x, order):
    return sum((F(v) * factorial(i) / factorial(i-order) * x**(i-order)
                for i, v in enumerate(coeff) if i >= order), F(0))

def main():
    if OUT.exists():
        raise FileExistsError(OUT)
    config = json.loads(CFG.read_text())
    fixtures_config = ROOT / config['fixtures_config']
    inputs = json.loads(fixtures_config.read_text())
    prior = json.loads((BASE / 'poly_split_diagnostic/summary.json').read_text())
    poly_path = ROOT / 'scripts/74_polynomial_forcing_prototype.py'
    oracle_path = ROOT / 'scripts/107_polynomial_derivative_crosscheck_20261101.py'
    if sha(poly_path) != prior['polynomial_source_sha256'] or sha(oracle_path) != prior['oracle_source_sha256'] or sha(fixtures_config) != prior['fixture_config_sha256']:
        raise ValueError('source/input hash mismatch')
    code = poly_path.read_text().split("fixture=ROOT/'configs/")[0]
    needle = "OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
    if code.count(needle) != 1:
        raise ValueError('unexpected prototype')
    ns = {'__file__': str(poly_path)}
    exec(compile(code.replace(needle, 'OUT=None'), str(poly_path), 'exec'), ns)
    Poly = ns['PolynomialForcing']
    oracle = load(oracle_path)
    class Omit(Poly):
        def source_derivative(self, x, j):
            return arb(0)
    class Reverse(Poly):
        def source_derivative(self, x, j):
            return -super().source_derivative(x, j)
    rows, blind = [], []
    points = [F(x) for x in config['points']]
    extra = F(config['extra_point'])
    for fixture in inputs['fixtures']:
        path = ROOT / fixture['path']
        if sha(path) != fixture['sha256']:
            raise ValueError('fixture digest mismatch')
        for net in json.loads(path.read_text())['networks']:
            args = (path, net['label'], fixture['sha256'], net['source_coefficients'], config['precision_bits'])
            problems = {'normal': Poly(*args), 'omit_source': Omit(*args), 'reverse_source_sign': Reverse(*args)}
            def value(p, x, order):
                fun = {0: p.residual, 2: p.r_second, 3: p.r_third, 4: p.r_fourth}[order]
                ball = fun(point(x))
                return exact_fraction(ball.lower()), exact_fraction(ball.upper())
            for x in points:
                for order in config['orders']:
                    r80, r110 = (F(*to_rational(oracle.mp_value(net, x, order, dps)._mpf_)) for dps in config['reference_digits'])
                    stable = abs(r80-r110) <= F(config['relative_stability']) * max(F(1), abs(r110))
                    src = source(net['source_coefficients'], x, order)
                    for fault, p in problems.items():
                        lo, hi = value(p, x, order)
                        rows.append({'fixture_sha256': fixture['sha256'], 'network': net['label'], 'source': net['source_name'], 'x': str(x), 'order': order, 'kind': fault, 'source_derivative': str(src), 'lower': str(lo), 'upper': str(hi), 'reference80': str(r80), 'reference110': str(r110), 'stable': stable, 'contained': lo <= r110 <= hi})
            # A smooth, explicitly injected residual error vanishing at all five
            # test points: detect it at an additional point, without changing Poly.
            for x in points + [extra]:
                r80, r110 = (F(*to_rational(oracle.mp_value(net, x, 0, dps)._mpf_)) for dps in config['reference_digits'])
                delta = F(1)
                for root in points:
                    delta *= x-root
                lo, hi = value(problems['normal'], x, 0)
                blind.append({'fixture_sha256': fixture['sha256'], 'network': net['label'], 'x': str(x), 'delta': str(delta), 'lower': str(lo+delta), 'upper': str(hi+delta), 'reference80': str(r80), 'reference110': str(r110), 'contained': lo+delta <= r110 <= hi+delta})
            print(net['label'], 'done', flush=True)
    BASE.mkdir(parents=True, exist_ok=True)
    stage_parent = Path(tempfile.mkdtemp(prefix='.multiorder_', dir=BASE))
    stage = stage_parent / 'bundle'
    stage.mkdir()
    (stage / 'locked_config.json').write_bytes(CFG.read_bytes())
    (stage / 'records.json').write_text(json.dumps(rows, indent=2)+'\n')
    (stage / 'blindspot.json').write_text(json.dumps(blind, indent=2)+'\n')
    summary = {'records': len(rows), 'blindspot_records': len(blind), 'stable': sum(r['stable'] for r in rows), 'normal_contained': sum(r['contained'] for r in rows if r['kind'] == 'normal'), 'nonzero_fault_detected': sum(not r['contained'] for r in rows if r['kind'] != 'normal' and F(r['source_derivative']) != 0), 'zero_fault_contained': sum(r['contained'] for r in rows if r['kind'] != 'normal' and F(r['source_derivative']) == 0), 'sampled_blindspot_contained': sum(r['contained'] for r in blind if F(r['x']) in points), 'extra_blindspot_detected': sum(not r['contained'] for r in blind if F(r['x']) == extra), 'records_sha256': sha(stage / 'records.json'), 'blindspot_sha256': sha(stage / 'blindspot.json'), 'config_sha256': sha(CFG), 'fixture_config_sha256': sha(fixtures_config), 'poly_sha256': sha(poly_path), 'oracle_sha256': sha(oracle_path), 'script_sha256': sha(Path(__file__)), 'scope': config['scope']}
    (stage / 'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    os.rename(stage, OUT)
    stage_parent.rmdir()
    print(json.dumps(summary), flush=True)
    return int(summary['records'] != 1440 or summary['blindspot_records'] != 144 or summary['stable'] != 1440 or summary['normal_contained'] != 480 or summary['nonzero_fault_detected'] != 416 or summary['zero_fault_contained'] != 544 or summary['sampled_blindspot_contained'] != 120 or summary['extra_blindspot_detected'] != 24)

if __name__ == '__main__':
    raise SystemExit(main())
