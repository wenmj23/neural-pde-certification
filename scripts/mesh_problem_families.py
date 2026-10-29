"""Manufactured off-center/multiple smooth peaks and independent references.

These are controlled error functions, not trained PINNs. Their maxima are not
assumed to equal the amplitude: a separate mpmath interval reference is used.
"""
from fractions import Fraction as F
import importlib.util
from pathlib import Path
from flint import arb
import mpmath as mp

from arb_cubic_reconstruction import BumpSix
from arb_problem_families import SmoothProblem
from arb_certification import finite, square, int_power
from interval_certification import context

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('mesh_reference18', ROOT/'scripts/18_continuous_reference_crosscheck.py')
ref18 = importlib.util.module_from_spec(spec); spec.loader.exec_module(ref18)


class ShiftedBump(BumpSix):
    def __init__(self, center, width, bits=103):
        self.center = F(center)
        if not 0 < self.center < 1 or not 0 < F(width) <= F(1, 2): raise ValueError('Invalid Gaussian parameters')
        super().__init__(width, bits)
        del self.known_error
        self.center_arb = arb(self.center.numerator)/self.center.denominator
        self.name = f'shifted_{center}_{width}'

    def _gs(self, x, order):
        t = (x-self.center_arb)/self.sigma; g = (-square(t)).exp(); values = []
        for j in range(order+1):
            p = arb(0)
            for c in reversed(self.coefficients[j]): p = p*t+c
            values.append(p*g/int_power(self.sigma, j))
        return values

    def mp_function(self):
        width = mp.mpf(self.width.numerator)/self.width.denominator
        center = mp.mpf(self.center.numerator)/self.center.denominator
        return lambda x: 4*x*(1-x)*mp.exp(-((x-center)/width)**2)/10000

    def error(self, x):
        self.guard()
        return finite(4*self.A*self.q(x)*(-square((x-self.center_arb)/self.sigma)).exp(), 'shifted manufactured error')


class GaussianMixture(SmoothProblem):
    def __init__(self, terms, bits=103):
        self.setup(bits)
        if not terms or any(F(t['weight']) < 0 for t in terms) or sum(F(t['weight']) for t in terms) != 1:
            raise ValueError('Nonnegative normalized Gaussian weights required')
        self.parts = [(F(t['weight']), ShiftedBump(t['center'], t['width'], bits)) for t in terms]

    def _evaluate(self, name, x):
        self.guard(); total = arb(0)
        for weight, part in self.parts:
            total += getattr(part, name)(x)*weight.numerator/weight.denominator
        return finite(total, 'manufactured mixture '+name)

    def residual(self, x): return self._evaluate('residual', x)
    def r_third(self, x): return self._evaluate('r_third', x)
    def r_fourth(self, x): return self._evaluate('r_fourth', x)
    def error(self, x): return self._evaluate('error', x)
    def mp_function(self):
        parts = [(mp.mpf(w.numerator)/w.denominator, p.mp_function()) for w, p in self.parts]
        return lambda x: sum(w*f(x) for w, f in parts)


class MPGaussianReference(ref18.MPReference):
    def __init__(self, terms):
        self.ctx = context(50); c = self.ctx
        def number(x):
            q = F(x); return c.mpf(q.numerator)/q.denominator
        self.terms = [(number(t['weight']), number(t['center']), number(t['width'])) for t in terms]
        self.node_cache = {}; self.residual_calls = 0

    def value_derivatives(self, x):
        c = self.ctx; value = first = second = c.mpf(0)
        for weight, center, width in self.terms:
            t = (x-center)/width; g = c.exp(-(t**2))
            value += weight*g
            first += weight*(-2*t/width)*g
            second += weight*((4*t**2-2)/width**2)*g
        return value, first, second

    def error(self, x):
        if x not in self.node_cache:
            xx = self.point(x); g, _, _ = self.value_derivatives(xx)
            self.node_cache[x] = 4*xx*(1-xx)*g/10000
        return self.node_cache[x]

    def residual(self, left, right):
        self.residual_calls += 1; c = self.ctx
        x = c.mpf([self.point(left).a, self.point(right).b])
        g, first, second = self.value_derivatives(x)
        return -4*(x*(1-x)*second+2*(1-2*x)*first-2*g)/10000
