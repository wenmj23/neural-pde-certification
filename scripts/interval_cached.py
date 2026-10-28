"""Coefficient caching for the reviewed fixed-weight interval implementation.

The coefficients a*w**j use exactly the same multiplication order and interval
context as before; only reuse changes. Precision cannot change after creation.
The conditional mathematical/interval-backend guarantees remain those of the
base class, and do not apply to arbitrary floating inference programs.
"""
from interval_certification import FrozenNeuralPoisson


class CachedFrozenNeuralPoisson(FrozenNeuralPoisson):
    def __init__(self, case_directory, ctx):
        super().__init__(case_directory, ctx)
        self._precision_bits = ctx.prec
        self._coefficients = [[a*w**j for j in range(5)] for w, a in zip(self.w, self.a)]

    def derivatives(self, x, order):
        if self.ctx.prec != self._precision_bits:
            raise ValueError('Rebuild cached coefficients after changing interval precision')
        if not 0 <= order <= 4:
            raise ValueError('Unsupported derivative order')
        ctx = self.ctx
        values = [self.c] + [ctx.mpf(0) for _ in range(order)]
        for w, b, coefficients in zip(self.w, self.b, self._coefficients):
            t = 1 - 2/(ctx.exp(2*(w*x+b))+1)
            ds = [t]
            if order >= 1: ds.append(1-t**2)
            if order >= 2: ds.append(-2*t+2*t**3)
            if order >= 3: ds.append(-2+8*t**2-6*t**4)
            if order >= 4: ds.append(16*t-40*t**3+24*t**5)
            for j in range(order+1):
                values[j] += coefficients[j]*ds[j]
        return values
