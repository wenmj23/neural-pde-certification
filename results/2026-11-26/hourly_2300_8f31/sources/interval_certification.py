"""Scoped interval certificate for fixed 1D Poisson candidates.

Conditional on correctness of mpmath 1.3.0's basic interval operations, exp,
sin and pi. The mathematical candidate has exact binary64 weights. This does
not certify arbitrary floating inference programs or general neural/PDE models.
See notes/working/interval_audit_2026-10-28.md for the proof and trust boundary.
"""
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import time

import mpmath
from mpmath import MPIntervalContext
from mpmath.libmp import to_rational

SUPPORTED_BACKEND = '1.3.0'


def context(dps=30):
    if mpmath.__version__ != SUPPORTED_BACKEND:
        raise RuntimeError('This audit targets mpmath==1.3.0; review other versions before use')
    if type(dps) is not int or not 25 <= dps <= 200:
        raise ValueError('Audit precision must be an integer from 25 to 200 decimal digits')
    ctx = MPIntervalContext()
    ctx.dps = dps
    return ctx


def endpoint_fraction(x, side=0):
    """Exact finite binary endpoint -> Fraction. No float/string conversion.

    _mpi_ and libmp.to_rational are isolated backend-specific interfaces; their
    semantics are tested and the backend version is pinned above.
    """
    endpoint = x._mpi_[side]
    if endpoint[3] < 0:
        raise ValueError('Nonfinite endpoint cannot be certified or exported')
    return Fraction(*to_rational(endpoint))


def export_endpoint(point, direction, digits=25):
    exact = endpoint_fraction(point)
    if exact != endpoint_fraction(point, 1):
        raise ValueError('Expected a point endpoint')
    if direction not in ('down', 'up'):
        raise ValueError('Unknown rounding direction')
    with localcontext() as dc:
        dc.prec = digits
        dc.rounding = ROUND_FLOOR if direction == 'down' else ROUND_CEILING
        text = str(Decimal(exact.numerator) / Decimal(exact.denominator))
    displayed = Fraction(text)
    if (direction == 'down' and displayed > exact) or (direction == 'up' and displayed < exact):
        raise ArithmeticError('Decimal export lost outward rounding')
    return dict(numerator=str(exact.numerator), denominator=str(exact.denominator),
                decimal_outward=text, decimal_rounding=direction)


def export_bounds(result):
    return dict(lower=export_endpoint(result['lower'], 'down'),
                upper=export_endpoint(result['upper'], 'up'))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FrozenNeuralPoisson:
    boundary_error = (0, 0)

    def __init__(self, case_directory, ctx):
        self.ctx = ctx
        directory = Path(case_directory)
        self.case = json.loads((directory / 'case.json').read_text())
        c = self.case
        n, p = c['network'], c['problem']
        if (c['schema_version'] != 1 or n['architecture'] != [1, 16, 1]
                or n['activation'] != 'tanh' or n['parameter_dtype'] != 'float64'
                or n['output_transform'] != 'x*(1-x)*raw_network(x)'
                or p['domain'] != [0., 1.] or p['operator'] != '-u_dd'
                or p['frequency'] != 1 or p['source'] != 'pi^2*sin(pi*x)'
                or p['boundary'] != {'type': 'Dirichlet', 'values': [0., 0.]}):
            raise ValueError('Unsupported candidate, operator, source or boundary conditions')
        for filename, digest in [(n['weights_file'], n['weights_sha256']),
                                 (c['provenance_file'], c['provenance_sha256']),
                                 (c['baseline_file'], c['baseline_sha256'])]:
            if sha256(directory / filename) != digest:
                raise ValueError('Frozen case file has changed: ' + filename)
        raw = json.loads((directory / n['weights_file']).read_text())
        expected_shapes = {'hidden.weight': (16, 1), 'hidden.bias': (16,),
                           'output.weight': (1, 16), 'output.bias': (1,)}
        if set(raw) != set(expected_shapes):
            raise ValueError('Unexpected weight names')
        for key, shape in expected_shapes.items():
            values = raw[key]
            if len(values) != shape[0] or (len(shape) == 2 and any(len(v) != shape[1] for v in values)):
                raise ValueError('Unexpected weight shape: ' + key)
        self.input_conversions = 0
        def exact_binary(value):
            if not isinstance(value, (float, int)) or isinstance(value, bool) or not math.isfinite(value):
                raise ValueError('Weight must be a finite binary64 number')
            binary = float(value)
            result = ctx.mpf(binary)
            target = Fraction.from_float(binary)
            if not endpoint_fraction(result) == endpoint_fraction(result, 1) == target:
                raise ArithmeticError('Weight was not imported as the exact frozen binary64 value')
            self.input_conversions += 1
            return result
        self.w = [exact_binary(v[0]) for v in raw['hidden.weight']]
        self.b = [exact_binary(v) for v in raw['hidden.bias']]
        self.a = [exact_binary(v) for v in raw['output.weight'][0]]
        self.c = exact_binary(raw['output.bias'][0])
        if self.input_conversions != n['parameter_count'] or self.input_conversions != 49:
            raise ValueError('Wrong parameter count')
        self.weights_sha256 = n['weights_sha256']
        self.omega = ctx.pi

    def derivatives(self, x, order):
        if not 0 <= order <= 4:
            raise ValueError('Unsupported derivative order')
        ctx = self.ctx
        values = [self.c] + [ctx.mpf(0) for _ in range(order)]
        for w, b, a in zip(self.w, self.b, self.a):
            # tanh has no direct iv implementation in the reviewed backend.
            t = 1 - 2 / (ctx.exp(2 * (w * x + b)) + 1)
            ds = [t]
            if order >= 1: ds.append(1 - t**2)
            if order >= 2: ds.append(-2*t + 2*t**3)
            if order >= 3: ds.append(-2 + 8*t**2 - 6*t**4)
            if order >= 4: ds.append(16*t - 40*t**3 + 24*t**5)
            for j in range(order + 1):
                values[j] += a * w**j * ds[j]
        return values

    def solution(self, x):
        return x * (1 - x) * self.derivatives(x, 0)[0]

    def residual(self, x):
        n, n1, n2 = self.derivatives(x, 2)
        return -x*(1-x)*n2 - 2*(1-2*x)*n1 + 2*n - self.omega**2*self.ctx.sin(self.omega*x)

    def r_second(self, x):
        _, _, n2, n3, n4 = self.derivatives(x, 4)
        return -x*(1-x)*n4 - 4*(1-2*x)*n3 + 12*n2 + self.omega**4*self.ctx.sin(self.omega*x)


def mesh(ctx, intervals):
    if type(intervals) is not int or intervals < 1:
        raise ValueError('Interval count must be a positive integer')
    nodes = [ctx.mpf(i) / intervals for i in range(intervals + 1)]
    if any(endpoint_fraction(b-a) <= 0 for a, b in zip(nodes[:-1], nodes[1:])):
        raise ValueError('Precision is insufficient to separate mesh nodes')
    return nodes


def signed_polynomial_norm(ctx, nodes, r):
    """Bound continuous cubic pieces through their Bernstein coefficient hull."""
    slope = ctx.mpf(0)
    for a, b, r0, r1 in zip(nodes[:-1], nodes[1:], r[:-1], r[1:]):
        h = b-a
        slope += h*(1-a)*(r0+r1)/2 - h**2*(r0+2*r1)/6
    value = ctx.mpf(0)
    lower, upper = ctx.mpf(0), ctx.mpf(0)
    for a, b, r0, r1 in zip(nodes[:-1], nodes[1:], r[:-1], r[1:]):
        h = b-a
        coefficients = [value, value+slope*h/3,
                        value+2*slope*h/3-r0*h**2/6,
                        value+slope*h-h**2*(2*r0+r1)/6]
        lower = max(lower, abs(coefficients[0]).a, abs(coefficients[-1]).a)
        upper = max(upper, *[abs(c).b for c in coefficients])
        value = coefficients[-1]
        slope -= h*(r0+r1)/2
    if 0 not in value:
        raise ArithmeticError('Auxiliary boundary endpoint failed its inclusion check')
    return lower, upper


def remainder_norm_upper(ctx, nodes, delta):
    slope = ctx.mpf(0)
    for a, b, d in zip(nodes[:-1], nodes[1:], delta):
        slope += d*(b-a)*(1-(a+b)/2)
    value, upper = ctx.mpf(0), ctx.mpf(0)
    for a, b, d in zip(nodes[:-1], nodes[1:], delta):
        h = b-a
        coefficients = [value, value+slope*h/2, value+slope*h-d*h**2/2]
        upper = max(upper, *[abs(c).b for c in coefficients])
        value = coefficients[-1]
        slope -= d*h
    if 0 not in value:
        raise ArithmeticError('Remainder boundary endpoint failed its inclusion check')
    return upper


def certify(problem, intervals):
    if problem.boundary_error != (0, 0):
        raise ValueError('This certificate requires exactly zero error at the boundaries')
    start = time.perf_counter()
    ctx = problem.ctx
    nodes = mesh(ctx, intervals)
    r = [problem.residual(x) for x in nodes]
    delta = []
    for a, b in zip(nodes[:-1], nodes[1:]):
        cell = ctx.mpf([a.a, b.b])
        K = abs(problem.r_second(cell)).b
        delta.append((K*(b-a)**2/8).b)
    zl, zu = signed_polynomial_norm(ctx, nodes, r)
    du = remainder_norm_upper(ctx, nodes, delta)
    result = dict(lower=max(ctx.mpf(0), (zl-du).a), upper=(zu+du).b,
                  z_lower=zl, z_upper=zu, remainder_upper=du,
                  intervals=intervals, precision_digits=ctx.dps,
                  elapsed_s=time.perf_counter()-start)
    # Export also rejects any nonfinite endpoints or reversed final range.
    if endpoint_fraction(result['lower']) > endpoint_fraction(result['upper']):
        raise ArithmeticError('Final bounds are reversed')
    export_bounds(result)
    return result


def continuous_reference(problem, intervals=2048):
    """Separate continuous error enclosure using known sin(pi*x), for audits.

    On each cell: |e| <= max(endpoint |e|) + sup|e''|*h^2/8, with e''=-r.
    This does not use r'' or the reconstructed auxiliary solution. The exact
    solution is used only here, never by certify(). Same interval backend.
    """
    start = time.perf_counter()
    ctx = problem.ctx
    nodes = mesh(ctx, intervals)
    errors = [problem.solution(x)-ctx.sin(ctx.pi*x) for x in nodes]
    lower = max(abs(e).a for e in errors)
    upper = ctx.mpf(0)
    for a, b, e0, e1 in zip(nodes[:-1], nodes[1:], errors[:-1], errors[1:]):
        cell = ctx.mpf([a.a, b.b])
        curvature = abs(problem.residual(cell)).b
        candidate = max(abs(e0).b, abs(e1).b) + curvature*(b-a)**2/8
        upper = max(upper, candidate.b)
    return dict(lower=lower, upper=upper, intervals=intervals,
                precision_digits=ctx.dps, elapsed_s=time.perf_counter()-start,
                method='endpoint error enclosures + per-cell interpolation remainder',
                scope='continuous error enclosure for the audit; shares the interval backend')
