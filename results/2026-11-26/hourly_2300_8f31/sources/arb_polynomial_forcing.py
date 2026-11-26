"""Research adapter for exact rational polynomial forcing, including nonintegers.

Public sine CLI remains unchanged. JSON floats denote their exact binary64
values; decimal/fraction strings denote exact rationals. Bounds depend on Arb.
"""
from fractions import Fraction
from math import factorial, isfinite
from flint import arb
from arb_cubic_reconstruction import SixthOrderDeep
from arb_certification import finite


def coefficient(value):
    if isinstance(value, bool):
        raise ValueError('boolean is not a source coefficient')
    if isinstance(value, Fraction):
        return value
    if type(value) is int:
        return Fraction(value)
    if type(value) is float:
        if not isfinite(value):
            raise ValueError('nonfinite source coefficient')
        return Fraction.from_float(value)
    if isinstance(value, str) and 0 < len(value) <= 120:
        try:
            return Fraction(value)
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError('invalid rational source coefficient') from exc
    raise ValueError('unsupported source coefficient')


class PolynomialForcing(SixthOrderDeep):
    def __init__(self, fixture, label, sha, coeff, bits=103):
        if not isinstance(coeff, (list, tuple)) or not coeff:
            raise ValueError('nonempty coefficient sequence required')
        self.coeff = tuple(coefficient(c) for c in coeff)
        super().__init__(fixture, label, sha, bits)

    def source_derivative(self, x, order):
        self.guard()
        if type(order) is not int or order < 0:
            raise ValueError('nonnegative integer derivative order required')
        out = arb(0)
        for i in reversed(range(order, len(self.coeff))):
            c = self.coeff[i] * (factorial(i) // factorial(i-order))
            out = out*x + arb(c.numerator)/c.denominator
        return finite(out, 'rational polynomial source derivative')

    def residual(self, x):
        n = self.network_jet(x, 2)
        return finite(-2*(self.q(x)*n[2]+(1-2*x)*n[1]-n[0])-self.source_derivative(x, 0), 'poly residual')

    def r_second(self, x):
        n = self.network_jet(x, 4)
        return finite(-24*(self.q(x)*n[4]+(1-2*x)*n[3]-n[2])-self.source_derivative(x, 2), 'poly residual second')

    def r_third(self, x):
        n = self.network_jet(x, 5)
        return finite(-120*(self.q(x)*n[5]+(1-2*x)*n[4]-n[3])-self.source_derivative(x, 3), 'poly residual third')

    def r_fourth(self, x):
        n = self.network_jet(x, 6)
        return finite(-720*(self.q(x)*n[6]+(1-2*x)*n[5]-n[4])-self.source_derivative(x, 4), 'poly residual fourth')
