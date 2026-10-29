"""Scoped smooth problem adapters for linear/quadratic residual experiments.

All candidates have zero error at x=0,1. The frozen sine-Poisson networks are
imported directly from a checksummed fixture; actual forcing frequency is
explicit. Analytic error families are controlled tests, not trained PINNs.
"""
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import flint
from flint import arb,ctx
import mpmath as mp
from arb_certification import finite,square,int_power,exact_fraction


def derivative_polys(t,order):
    ds=[t]
    if order==0:return ds
    t2=square(t);s=1-t2;ds.append(s)
    if order>=2:ds.append(-2*t*s)
    if order>=3:ds.append(-2*s*(1-3*t2))
    if order>=4:ds.append(8*t*s*(2-3*t2))
    if order>=5:ds.append(8*s*(2+t2*(-15+15*t2)))
    return ds


class SmoothProblem:
    boundary_error=(0,0)
    def setup(self,bits):
        if flint.__version__!='0.9.0':raise RuntimeError('Review backend version before use')
        if type(bits) is not int or not 80<=bits<=700:raise ValueError('Unsupported precision')
        ctx.prec=bits;ctx.threads=1;self.bits=bits
    def guard(self):
        if ctx.prec!=self.bits:raise ValueError('Precision changed; rebuild problem/cache')
    def q(self,x):return arb(1)/4-square(x-arb(1)/2)


class FrozenSinePoisson(SmoothProblem):
    def __init__(self,fixture,label,expected_sha256,bits=103):
        self.setup(bits);data=Path(fixture).read_bytes()
        if hashlib.sha256(data).hexdigest()!=expected_sha256:raise ValueError('Frozen fixture changed')
        rows=json.loads(data)['networks'];matches=[r for r in rows if r['label']==label]
        if len(matches)!=1:raise ValueError('Unknown/ambiguous frozen network')
        row=matches[0];k=row['frequency']
        if type(k) is not int or k<1:raise ValueError('Positive integer sine frequency required')
        self.name=label;self.frequency=k;self.seed=row['seed'];self.raw=row['weights']
        shapes={'hidden.weight':(16,1),'hidden.bias':(16,),'output.weight':(1,16),'output.bias':(1,)}
        if set(self.raw)!=set(shapes):raise ValueError('Unexpected weight names')
        for key,shape in shapes.items():
            values=self.raw[key]
            if len(values)!=shape[0] or (len(shape)==2 and any(len(v)!=shape[1] for v in values)):
                raise ValueError('Wrong weight shape')
        def value(v):
            if isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v):raise ValueError('Invalid weight')
            n,d=float(v).as_integer_ratio();result=arb(n)/d
            if exact_fraction(result)!=Fraction(n,d):raise ArithmeticError('Inexact weight import')
            return result
        self.w=[value(v[0]) for v in self.raw['hidden.weight']];self.b=[value(v) for v in self.raw['hidden.bias']]
        self.a=[value(v) for v in self.raw['output.weight'][0]];self.c=value(self.raw['output.bias'][0])
        self.coeff=[[a*int_power(w,j) for j in range(6)] for a,w in zip(self.a,self.w)]
        self.omega=k*arb.pi();self.wp=[int_power(self.omega,j) for j in range(6)]
        self.identity=dict(fixture_sha256=expected_sha256,label=label,forcing='(k*pi)^2*sin(k*pi*x)',frequency=k,
            weights_sha256=hashlib.sha256(json.dumps(self.raw,sort_keys=True).encode()).hexdigest(),
            historical_grid_error=row['reference_error_grid'],reference_kind='observed grid error only')
    def derivatives(self,x,orders):
        self.guard();result={j:(self.c if j==0 else arb(0)) for j in orders}
        for w,b,co in zip(self.w,self.b,self.coeff):
            t=finite((w*x+b).tanh(),'tanh');ds=derivative_polys(t,max(orders))
            for j in orders:result[j]+=co[j]*ds[j]
        return result
    def solution(self,x):return self.q(x)*self.derivatives(x,(0,))[0]
    def residual(self,x):
        d=self.derivatives(x,(0,1,2))
        return -self.q(x)*d[2]-2*(1-2*x)*d[1]+2*d[0]-self.wp[2]*(self.omega*x).sin()
    def r_second(self,x):
        d=self.derivatives(x,(2,3,4))
        return -self.q(x)*d[4]-4*(1-2*x)*d[3]+12*d[2]+self.wp[4]*(self.omega*x).sin()
    def r_third(self,x):
        d=self.derivatives(x,(3,4,5))
        return -self.q(x)*d[5]-5*(1-2*x)*d[4]+20*d[3]+self.wp[5]*(self.omega*x).cos()
    def error(self,x):return self.solution(x)-(self.omega*x).sin()
    def mp_function(self):
        w=[mp.mpf(float(v[0])) for v in self.raw['hidden.weight']];b=[mp.mpf(float(v)) for v in self.raw['hidden.bias']]
        a=[mp.mpf(float(v)) for v in self.raw['output.weight'][0]];c=mp.mpf(float(self.raw['output.bias'][0]));k=self.frequency
        return lambda x:x*(1-x)*(c+sum(aa*mp.tanh(ww*x+bb) for ww,bb,aa in zip(w,b,a)))-mp.sin(k*mp.pi*x)


class SineError(SmoothProblem):
    def __init__(self,frequency,bits=103,frequencies=None,coefficients=None,name=None):
        self.setup(bits);self.frequencies=frequencies or [frequency];self.coefficients=coefficients or [Fraction(1)]
        if len(self.frequencies)!=len(self.coefficients) or any(type(k) is not int or k<1 for k in self.frequencies):
            raise ValueError('Matching coefficients and positive integer frequencies required')
        if sum(abs(c) for c in self.coefficients)!=1:raise ValueError('Coefficients must normalize known maximum')
        if any(k%2!=1 for k in self.frequencies) and len(self.frequencies)>1:raise ValueError('Mixed test needs odd frequencies')
        if len(self.frequencies)>1 and any(c*(-1)**((k-1)//2)<0 for k,c in zip(self.frequencies,self.coefficients)):
            raise ValueError('Mixed components must align at midpoint')
        self.name=name or f'sine_{frequency}';self.known_error=Fraction(1,10000);self.A=arb(1)/10000
        self.omegas=[k*arb.pi() for k in self.frequencies]
        self.co=[arb(c.numerator)/c.denominator for c in self.coefficients]
    def _r(self,x,order):
        self.guard();result=arb(0)
        for w,c in zip(self.omegas,self.co):
            trig=(w*x).cos() if order%2 else (w*x).sin()
            sign=1 if order==0 else -1
            result+=sign*self.A*c*int_power(w,order+2)*trig
        return result
    def residual(self,x):return self._r(x,0)
    def r_second(self,x):return self._r(x,2)
    def r_third(self,x):return self._r(x,3)
    def error(self,x):return self.A*sum(c*(w*x).sin() for w,c in zip(self.omegas,self.co))
    def mp_function(self):
        terms=[(k,mp.mpf(c.numerator)/c.denominator) for k,c in zip(self.frequencies,self.coefficients)]
        return lambda x:sum(c*mp.sin(k*mp.pi*x) for k,c in terms)/10000


class BumpError(SmoothProblem):
    def __init__(self,width,bits=103):
        self.setup(bits);self.width=Fraction(str(width));self.sigma=arb(self.width.numerator)/self.width.denominator
        self.A=arb(1)/10000;self.known_error=Fraction(1,10000);self.name=f'bump_{width}'
        self.coefficients=[[1],[0,-2],[-2,0,4],[0,12,0,-8],[12,0,-48,0,16],[0,-120,0,160,0,-32]]
    def _gs(self,x,order):
        t=(x-arb(1)/2)/self.sigma;g=(-square(t)).exp();values=[]
        for j in range(order+1):
            p=arb(0)
            for c in reversed(self.coefficients[j]):p=p*t+c
            values.append(p*g/int_power(self.sigma,j))
        return values
    def _e_derivative(self,x,order):
        self.guard();g=self._gs(x,order)
        return 4*self.A*(self.q(x)*g[order]+order*(1-2*x)*g[order-1]-order*(order-1)*g[order-2])
    def residual(self,x):return -self._e_derivative(x,2)
    def r_second(self,x):return -self._e_derivative(x,4)
    def r_third(self,x):return -self._e_derivative(x,5)
    def error(self,x):return 4*self.A*self.q(x)*(-square((x-arb(1)/2)/self.sigma)).exp()
    def mp_function(self):
        s=mp.mpf(self.width.numerator)/self.width.denominator
        return lambda x:4*x*(1-x)*mp.exp(-((x-mp.mpf(1)/2)/s)**2)/10000


class QuarticError(SmoothProblem):
    def __init__(self,bits=103):
        self.setup(bits);self.name='quartic_exact_quadratic_residual';self.A=arb(1)/10000;self.known_error=Fraction(1,10000)
    def residual(self,x):return self.A*(-32+192*x-192*square(x))
    def r_second(self,x):return -384*self.A
    def r_third(self,x):return arb(0)
    def error(self,x):return 16*self.A*square(self.q(x))
    def mp_function(self):return lambda x:16*x*x*(1-x)**2/10000
