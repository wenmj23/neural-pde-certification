"""Interval expression experiments, keeping the audited Poisson model fixed.

Variants: centered q, factored tanh derivative polynomials, and stationary-point
polynomial ranges. All operations are interval operations in the pinned backend.
No sampled derivative maxima are used as bounds. Each variant inherits the base
case/hash/precision guards and the same conditional guarantee scope.
"""
from interval_cached import CachedFrozenNeuralPoisson


def factored_polynomials(t,order=4):
    s=1-t**2
    ds=[t,s]
    if order>=2:ds.append(-2*t*s)
    if order>=3:ds.append(-2*s*(1-3*t**2))
    if order>=4:ds.append(8*t*s*(2-3*t**2))
    return ds


class ReformulatedPoisson(CachedFrozenNeuralPoisson):
    def __init__(self,case_directory,ctx,strategy='factored'):
        if strategy not in ('centered','factored','critical'):
            raise ValueError('Unsupported interval expression strategy')
        super().__init__(case_directory,ctx)
        self.strategy=strategy
        self._critical={}
        if strategy=='critical':
            s13=ctx.sqrt(ctx.mpf(1)/3)
            s23=ctx.sqrt(ctx.mpf(2)/3)
            low=ctx.sqrt((15-ctx.sqrt(105))/30)
            high=ctx.sqrt((15+ctx.sqrt(105))/30)
            roots={1:[ctx.mpf(0)],2:[-s13,s13],3:[-s23,ctx.mpf(0),s23],
                   4:[-high,-low,low,high]}
            for order,points in roots.items():
                self._critical[order]=[(p,factored_polynomials(p)[order]) for p in points]

    def _tanh(self,z):
        c=self.ctx
        t=1-2/(c.exp(2*z)+1)
        lo,hi=max(t.a,c.mpf(-1)),min(t.b,c.mpf(1))
        if lo>hi:raise ArithmeticError('Empty tanh range')
        return c.mpf([lo,hi])

    def _range_poly(self,t,order):
        values=[factored_polynomials(t.a,order)[order],factored_polynomials(t.b,order)[order]]
        for point,value in self._critical.get(order,[]):
            # Include any stationary point whose enclosure might meet [t.a,t.b].
            if point.b>=t.a and point.a<=t.b:values.append(value)
        return self.ctx.mpf([min(v.a for v in values),max(v.b for v in values)])

    def _selected(self,x,orders,tight=False):
        if self.ctx.prec!=self._precision_bits:
            raise ValueError('Rebuild after changing precision')
        c=self.ctx
        result={j:(self.c if j==0 else c.mpf(0)) for j in orders}
        for w,b,coeff in zip(self.w,self.b,self._coefficients):
            t=self._tanh(w*x+b)
            if tight and self.strategy=='critical':
                ds={j:self._range_poly(t,j) for j in orders}
            else:
                all_ds=factored_polynomials(t,max(orders))
                ds={j:all_ds[j] for j in orders}
            for j in orders:result[j]+=coeff[j]*ds[j]
        return result

    def derivatives(self,x,order):
        if self.strategy=='centered':return super().derivatives(x,order)
        if not 0<=order<=4:raise ValueError('Unsupported derivative order')
        values=self._selected(x,range(order+1))
        return [values[j] for j in range(order+1)]

    def q(self,x):return self.ctx.mpf(1)/4-(x-self.ctx.mpf(1)/2)**2

    def solution(self,x):return self.q(x)*self.derivatives(x,0)[0]

    def residual(self,x):
        n,n1,n2=self.derivatives(x,2)
        return -self.q(x)*n2-2*(1-2*x)*n1+2*n-self.omega**2*self.ctx.sin(self.omega*x)

    def r_second(self,x):
        if self.strategy=='centered':
            _,_,n2,n3,n4=self.derivatives(x,4)
        else:
            d=self._selected(x,(2,3,4),tight=True)
            n2,n3,n4=d[2],d[3],d[4]
        return -self.q(x)*n4-4*(1-2*x)*n3+12*n2+self.omega**4*self.ctx.sin(self.omega*x)
