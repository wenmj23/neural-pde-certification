"""New analytic error families for the existing 1D Poisson certificate.

The true solution is u=x(1-x), and the candidate is v=u+e. Thus r=-e''.
All nonzero examples have analytically known ||e||_infinity=|A|. The verifier
receives only residual callbacks and derivative bounds, never the true error.
Ordinary float64 experiments, not outward-rounded validated numerics.
"""
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import platform
import statistics
import sys
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpmath as mp
from certification_core import certificate, adaptive_nodes, trig_range

ROOT = Path(__file__).resolve().parents[1]
METHODS = [('constant', False), ('linear', False), ('linear', True)]
BUDGETS = [16, 32, 64, 128, 256, 512, 1000, 2000]


@dataclass
class MixedSine:
    name: str
    frequencies: tuple
    coefficients: tuple
    amplitude: float = 1e-4

    def __post_init__(self):
        assert len(self.frequencies) == len(self.coefficients)
        assert all(k > 0 and k % 2 == 1 for k in self.frequencies)
        assert sum(abs(c) for c in self.coefficients) == 1
        # Every component reaches its positive signed peak at x=1/2.
        assert all(c * (-1) ** ((k - 1) // 2) >= 0
                   for k, c in zip(self.frequencies, self.coefficients))

    @property
    def exact_error(self): return abs(self.amplitude)

    @property
    def peak(self): return .5

    def error(self, x):
        return self.amplitude * sum(c * np.sin(k * np.pi * x)
                                   for k, c in zip(self.frequencies, self.coefficients))

    def error_mp(self, x):
        return mp.mpf(self.amplitude) * sum(mp.mpf(c) * mp.sin(k * mp.pi * x)
                                          for k, c in zip(self.frequencies, self.coefficients))

    def residual(self, x):
        return self.amplitude * sum(c * (k * np.pi) ** 2 * np.sin(k * np.pi * x)
                                   for k, c in zip(self.frequencies, self.coefficients))

    def derivative_bound(self, nodes, order):
        lower = np.zeros(len(nodes) - 1)
        upper = lower.copy()
        for k, c in zip(self.frequencies, self.coefficients):
            lo, hi = trig_range(nodes[:-1], nodes[1:], k * np.pi,
                                'cos' if order == 1 else 'sin')
            coefficient = self.amplitude * c * (k * np.pi) ** (order + 2)
            if order == 2: coefficient = -coefficient
            lower += np.minimum(coefficient * lo, coefficient * hi)
            upper += np.maximum(coefficient * lo, coefficient * hi)
        return np.maximum(abs(lower), abs(upper))


def falling(n, j):
    return math.prod(range(n - j + 1, n + 1))


def beta_shape_range(left, right, a, b):
    """Exact-real extrema of x^a(1-x)^b on a cell, evaluated in float64."""
    v0 = left ** a * (1 - left) ** b
    v1 = right ** a * (1 - right) ** b
    lo, hi = np.minimum(v0, v1), np.maximum(v0, v1)
    if a + b:
        peak = a / (a + b)
        hi = np.where((left <= peak) & (peak <= right),
                      peak ** a * (1 - peak) ** b, hi)
    return lo, hi


@dataclass
class BetaError:
    name: str
    p: int
    q: int
    amplitude: float = 1e-4

    @property
    def peak(self): return self.p / (self.p + self.q)

    @property
    def normalization(self): return self.peak ** self.p * (1 - self.peak) ** self.q

    @property
    def exact_error(self): return abs(self.amplitude)

    def terms(self, derivative):
        for j in range(derivative + 1):
            if j > self.p or derivative - j > self.q: continue
            coefficient = (math.comb(derivative, j) * falling(self.p, j)
                           * falling(self.q, derivative - j) * (-1) ** (derivative - j))
            yield coefficient, self.p - j, self.q - derivative + j

    def error(self, x):
        return self.amplitude / self.normalization * x ** self.p * (1 - x) ** self.q

    def error_mp(self, x):
        peak = mp.mpf(self.p) / (self.p + self.q)
        norm = peak ** self.p * (1 - peak) ** self.q
        return mp.mpf(self.amplitude) / norm * x ** self.p * (1 - x) ** self.q

    def residual(self, x):
        return -self.amplitude / self.normalization * sum(
            c * x ** a * (1 - x) ** b for c, a, b in self.terms(2))

    def derivative_bound(self, nodes, order):
        lower = np.zeros(len(nodes) - 1)
        upper = lower.copy()
        for c, a, b in self.terms(order + 2):
            lo, hi = beta_shape_range(nodes[:-1], nodes[1:], a, b)
            c = -self.amplitude / self.normalization * c
            lower += np.minimum(c * lo, c * hi)
            upper += np.maximum(c * lo, c * hi)
        return np.maximum(abs(lower), abs(upper))


@dataclass
class AffineResidualError:
    name: str = 'affine_residual'
    amplitude: float = 1e-4

    @property
    def exact_error(self): return abs(self.amplitude)

    @property
    def peak(self): return 1 / np.sqrt(3)

    @property
    def scale(self): return self.amplitude * 3 * np.sqrt(3) / 2

    def error(self, x): return self.scale * (x - x ** 3)

    def error_mp(self, x):
        return mp.mpf(self.amplitude) * 3 * mp.sqrt(3) / 2 * (x - x ** 3)

    def residual(self, x): return 6 * self.scale * x

    def derivative_bound(self, nodes, order):
        return np.full(len(nodes) - 1, abs(6 * self.scale) if order == 1 else 0.)


def problems():
    return [
        AffineResidualError(),
        MixedSine('mixed_mild', (1, 5, 21), (.75, .1875, .0625)),
        MixedSine('mixed_signed', (1, 3, 41), (.5, -.25, .25)),
        MixedSine('weak_high_tail', (1, 5, 81), (.9375, .0546875, .0078125)),
        BetaError('broad_polynomial', 2, 2),
        BetaError('left_peak_1_16', 1, 15),
        BetaError('left_peak_1_64', 1, 63),
        BetaError('right_peak_63_64', 63, 1),
        MixedSine('zero_error', (1,), (1.,), 0.),
    ]


def check_adapters(cases):
    """Independent high-precision differentiation, used to catch formula errors."""
    mp.mp.dps = 60
    rng = np.random.default_rng(20261028)
    nodes = np.r_[0., np.sort(rng.uniform(size=72)), 1.]
    grid = np.linspace(0, 1, 10001)
    for problem in cases:
        amplitude = problem.exact_error
        assert np.max(abs(problem.error(grid))) <= amplitude * (1 + 2e-13) + 1e-30
        assert abs(abs(problem.error(problem.peak)) - amplitude) <= amplitude * 2e-13 + 1e-30
        assert np.max(abs(problem.error(np.array([0., 1.])))) <= amplitude * 1e-12 + 1e-30
        for xx in (.003, .017, .031, .13, .37, .5, .79, .983, .997):
            x = mp.mpf(xx)
            reference = float(-mp.diff(problem.error_mp, x, 2))
            value = float(problem.residual(np.array(xx)))
            assert np.isclose(value, reference, rtol=2e-10, atol=1e-12), (problem.name, xx)
            cell = min(np.searchsorted(nodes, xx, side='right') - 1, len(nodes) - 2)
            for order in (1, 2):
                derivative = abs(float(-mp.diff(problem.error_mp, x, order + 2)))
                bound = problem.derivative_bound(nodes, order)[cell]
                assert derivative <= bound * (1 + 2e-10) + 1e-11, (problem.name, xx, order)
    # A linear residual has identically zero second derivative; its reconstruction
    # is exact in real arithmetic on any mesh, including nonuniform meshes.
    result = certificate(cases[0], nodes, 'linear')
    assert result['D'] == 0
    assert np.isclose(result['upper'], cases[0].exact_error, rtol=2e-12, atol=0)
    print('Analytic peaks, boundary values, residual formulas, derivative enclosures and affine exactness checked.', flush=True)


def run_case(problem, n, method, adaptive):
    if adaptive:
        nodes, mesh = adaptive_nodes(problem, n, method)
    else:
        nodes = np.linspace(0, 1, n + 1)
        mesh = dict(mesh_ms=0., bound_evaluations=0)
    result = certificate(problem, nodes, method)
    exact = problem.exact_error
    # Relative float64 tolerance only; these assertions are numerical checks.
    slack = exact * 2e-9 + 1e-30
    assert result['lower'] <= exact + slack, (problem.name, n, method, adaptive, 'lower')
    assert result['upper'] >= exact - slack, (problem.name, n, method, adaptive, 'upper')
    if exact == 0: assert result['lower'] == result['upper'] == 0
    return dict(problem=problem.name, intervals=n, method=method, adaptive=adaptive,
                analytic_error=exact, lower=result['lower'], upper=result['upper'],
                ratio=result['upper'] / exact if exact else None, remainder=result['D'],
                total_ms=result['elapsed_ms'] + mesh['mesh_ms'],
                interval_bound_evaluations=n + mesh['bound_evaluations'])


def matched_benchmark(cases):
    spec = importlib.util.spec_from_file_location('matched08', Path(__file__).with_name('08_matched_tolerance_benchmark.py'))
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    output = []
    for problem in cases:
        if problem.exact_error == 0: continue
        for method, adaptive in METHODS:
            runs = [helper.solve(problem, method, adaptive, cap=65536) for _ in range(3)]
            r = runs[0]
            assert all(x['intervals'] == r['intervals'] and x['success'] == r['success'] for x in runs)
            r.update(problem=problem.name, method=method, adaptive=adaptive,
                     median_total_ms=statistics.median(x['total_ms'] for x in runs),
                     timing_samples_ms=[x['total_ms'] for x in runs],
                     ratio_to_analytic_error=r['upper'] / problem.exact_error)
            output.append(r)
    return output


def make_plot(cases, records, out):
    colors = ['#577590', '#43aa8b', '#f3722c']
    labels = ['Constant / uniform', 'Linear / uniform', 'Linear / adaptive']
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.8))
    selected = ['mixed_signed', 'weak_high_tail', 'left_peak_1_64']
    x = np.linspace(0, 1, 10001)
    for j, name in enumerate(selected):
        p = next(p for p in cases if p.name == name)
        axes[0, j].plot(x, p.error(x) / p.exact_error, color='#277da1')
        axes[0, j].set(title=name, xlabel='x', ylabel='Error / known maximum')
        for (method, adaptive), label, color in zip(METHODS, labels, colors):
            rows = [r for r in records if r['problem'] == name and r['method'] == method and r['adaptive'] == adaptive]
            axes[1, j].loglog([r['intervals'] for r in rows], [r['ratio'] for r in rows], '-o', ms=3, color=color, label=label)
        axes[1, j].axhline(1, color='gray', lw=.8, ls='--')
        axes[1, j].set(xlabel='Final intervals', ylabel='Upper bound / analytic error')
        axes[1, j].legend(fontsize=7)
    for ax in axes.flat: ax.grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(out / 'comparison.png', dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    out = args.out or ROOT / 'results/generated/10_function_generalization' / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    # Each run is a fresh evidence directory; never overwrite previous outputs.
    out.mkdir(parents=True, exist_ok=False)
    cases = problems()
    config = dict(pde='-u_dd=2 on [0,1], zero Dirichlet; u=x(1-x), v=u+e, r=-e_dd',
                  cases=[dict(type=type(p).__name__, **asdict(p)) for p in cases],
                  budgets=BUDGETS, methods=METHODS, matched_gap_target=1.2,
                  matched_cap=65536, matched_repeats=3, arithmetic='ordinary float64',
                  python=sys.version, platform=platform.platform(), numpy=np.__version__,
                  matplotlib=matplotlib.__version__, mpmath=mp.__version__)
    (out / 'config.json').write_text(json.dumps(config, indent=2) + '\n')
    check_adapters(cases)
    records = []
    for problem in cases:
        for n in BUDGETS:
            for method, adaptive in METHODS:
                records.append(run_case(problem, n, method, adaptive))
        print(problem.name, [(r['method'], r['adaptive'], round(r['ratio'], 8) if r['ratio'] is not None else None)
                              for r in records if r['problem'] == problem.name and r['intervals'] == 1000], flush=True)
    for n in BUDGETS:
        for method in ('constant', 'linear'):
            a = next(r for r in records if r['problem'] == 'left_peak_1_64' and r['intervals'] == n and r['method'] == method and not r['adaptive'])
            b = next(r for r in records if r['problem'] == 'right_peak_63_64' and r['intervals'] == n and r['method'] == method and not r['adaptive'])
            assert np.isclose(a['upper'], b['upper'], rtol=1e-7, atol=1e-14)
    (out / 'records.json').write_text(json.dumps(records, indent=2) + '\n')
    with (out / 'records.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    matched = matched_benchmark(cases)
    (out / 'matched_benchmark.json').write_text(json.dumps(matched, indent=2) + '\n')
    make_plot(cases, records, out)
    sources = [Path(__file__), Path(__file__).with_name('certification_core.py'),
               Path(__file__).with_name('08_matched_tolerance_benchmark.py')]
    (out / 'source_hashes.json').write_text(json.dumps(
        {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}, indent=2) + '\n')
    print('All checks passed:', len(records), 'fixed-budget configurations;', len(matched), 'matched-gap configurations.', flush=True)
    print('RESULT_DIRECTORY', out, flush=True)


if __name__ == '__main__': main()
