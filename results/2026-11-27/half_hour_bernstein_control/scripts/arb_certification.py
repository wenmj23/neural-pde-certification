"""Native FLINT/Arb backend for the same scoped 1D certificate formulas.

Uses exact binary64 weights, outward ball arithmetic, Bernstein hulls and exact
rational endpoint export. Requires python-flint==0.9.0. Global FLINT precision
is guarded; use sequentially, not concurrently with other precision settings.
This backend port has its own review/tests, not the status of the old audit.
"""
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
from fractions import Fraction
import json
from pathlib import Path
import time

import flint
from flint import arb,ctx
from interval_certification import FrozenNeuralPoisson,context as mp_context


def finite(x,label='interval arithmetic'):
    # Arb's indeterminate values do not behave like ordinary ordered reals.
    # Reject before max/min/abs_lower can silently discard an invalid result.
    if not x.is_finite():raise ArithmeticError('Nonfinite value in '+label)
    return x


def exact_fraction(point):
    if not point.is_exact() or not point.is_finite():
        raise ValueError('Need a finite exact endpoint')
    q=point.fmpq()
    return Fraction(int(q.p),int(q.q))


def export(point,direction,digits=25):
    q=exact_fraction(point)
    if direction not in ('down','up'):raise ValueError('Unknown rounding direction')
    with localcontext() as c:
        c.prec=digits;c.rounding=ROUND_FLOOR if direction=='down' else ROUND_CEILING
        text=str(Decimal(q.numerator)/Decimal(q.denominator))
    d=Fraction(text)
    if (direction=='down' and d>q) or (direction=='up' and d<q):
        raise ArithmeticError('Unsafe decimal rounding')
    return dict(numerator=str(q.numerator),denominator=str(q.denominator),
                decimal_outward=text,decimal_rounding=direction)


def bounds(result):
    return dict(lower=export(result['lower'],'down'),upper=export(result['upper'],'up'))


def square(x):
    # arb ** 2 returned indeterminate for zero-crossing balls in the tested
    # Python binding. Use abs endpoint bounds and ordinary multiplication.
    finite(x,'square input')
    a,b=x.abs_lower(),x.abs_upper()
    return finite((a*a).lower().union((b*b).upper()),'square output')


def polys(t,order=4):
    t2=square(t);s=1-t2
    ds=[t,s]
    if order>=2:ds.append(-2*t*s)
    if order>=3:ds.append(-2*s*(1-3*t2))
    if order>=4:ds.append(8*t*s*(2-3*t2))
    return ds


class ArbPoisson:
    boundary_error=(0,0)
    def __init__(self,case_directory,bits=103,strategy='factored'):
        if flint.__version__!='0.9.0':raise RuntimeError('Review a new python-flint version first')
        if type(bits) is not int or not 80<=bits<=700:raise ValueError('Unsupported precision')
        if strategy not in ('expanded','factored','critical'):raise ValueError('Unknown strategy')
        # Reuse the existing metadata/hash/shape guards; import into Arb afresh.
        reviewed=FrozenNeuralPoisson(case_directory,mp_context(30))
        self.case=reviewed.case;self.weights_sha256=reviewed.weights_sha256
        ctx.prec=bits;ctx.threads=1;self.bits=bits;self.strategy=strategy
        raw=json.loads((Path(case_directory)/self.case['network']['weights_file']).read_text())
        def value(v):
            numerator,denominator=float(v).as_integer_ratio()
            a=arb(numerator)/denominator
            if exact_fraction(a)!=Fraction(numerator,denominator):
                raise ArithmeticError('Frozen binary64 input not exact')
            return a
        self.w=[value(v[0]) for v in raw['hidden.weight']]
        self.b=[value(v) for v in raw['hidden.bias']]
        self.a=[value(v) for v in raw['output.weight'][0]]
        self.c=value(raw['output.bias'][0]);self.omega=arb.pi()
        self.coeff=[[a*int_power(w,j) for j in range(5)] for w,a in zip(self.w,self.a)]
        self.critical={}
        if strategy=='critical':
            s13=(arb(1)/3).sqrt();s23=(arb(2)/3).sqrt()
            low=((15-arb(105).sqrt())/30).sqrt();high=((15+arb(105).sqrt())/30).sqrt()
            roots={2:[-s13,s13],3:[-s23,arb(0),s23],4:[-high,-low,low,high]}
            self.critical={j:[(v,polys(v)[j]) for v in vs] for j,vs in roots.items()}

    def guard(self):
        if ctx.prec!=self.bits:raise ValueError('FLINT precision changed; rebuild candidate')

    def _range(self,t,j):
        values=[polys(t.lower(),j)[j],polys(t.upper(),j)[j]]
        for p,v in self.critical.get(j,[]):
            if p.upper()>=t.lower() and p.lower()<=t.upper():values.append(v)
        lo=min(v.lower() for v in values);hi=max(v.upper() for v in values)
        return lo.union(hi)

    def derivatives(self,x,orders,tight=False):
        self.guard();values={j:(self.c if j==0 else arb(0)) for j in orders}
        for w,b,coeff in zip(self.w,self.b,self.coeff):
            t=finite((w*x+b).tanh(),'tanh')
            if tight and self.strategy=='critical':
                ds={j:self._range(t,j) for j in orders}
            elif self.strategy=='expanded':
                t2=square(t);t3=t*t2;t4=square(t2);t5=t*t4
                all_ds=[t,1-t2,-2*t+2*t3,-2+8*t2-6*t4,16*t-40*t3+24*t5]
                ds={j:all_ds[j] for j in orders}
            else:
                all_ds=polys(t,max(orders));ds={j:all_ds[j] for j in orders}
            for j in orders:values[j]+=coeff[j]*ds[j]
        return values

    def q(self,x):
        return x*(1-x) if self.strategy=='expanded' else arb(1)/4-square(x-arb(1)/2)

    def solution(self,x):return self.q(x)*self.derivatives(x,(0,))[0]

    def residual(self,x):
        n=self.derivatives(x,(0,1,2))
        return -self.q(x)*n[2]-2*(1-2*x)*n[1]+2*n[0]-self.omega*self.omega*(self.omega*x).sin()

    def r_second(self,x):
        d=self.derivatives(x,(2,3,4),tight=True)
        return -self.q(x)*d[4]-4*(1-2*x)*d[3]+12*d[2]+int_power(self.omega,4)*(self.omega*x).sin()


def int_power(x,n):
    if n<0 or type(n) is not int:raise ValueError('Nonnegative integer power required')
    if n==0:return arb(1)
    if n==1:return x
    if n%2==0:return square(int_power(x,n//2))
    return x*int_power(x,n-1)


def make_mesh(n):
    if type(n) is not int or n<1:raise ValueError('Positive interval count required')
    nodes=[finite(arb(i)/n,'mesh node') for i in range(n+1)]
    if any((b-a).lower()<=0 for a,b in zip(nodes[:-1],nodes[1:])):
        raise ValueError('Unseparated mesh nodes')
    return nodes


def polynomial_norm(nodes,r):
    slope=arb(0)
    for a,b,r0,r1 in zip(nodes[:-1],nodes[1:],r[:-1],r[1:]):
        h=b-a
        slope=finite(slope+h*(1-a)*(r0+r1)/2-h*h*(r0+2*r1)/6,'signed initial slope')
    value=arb(0);lo=arb(0);hi=arb(0)
    for a,b,r0,r1 in zip(nodes[:-1],nodes[1:],r[:-1],r[1:]):
        h=b-a
        coeff=[value,value+slope*h/3,value+2*slope*h/3-r0*h*h/6,
               value+slope*h-h*h*(2*r0+r1)/6]
        for v in coeff:finite(v,'signed Bernstein coefficient')
        lo=max(lo,coeff[0].abs_lower(),coeff[-1].abs_lower())
        hi=max(hi,*[v.abs_upper() for v in coeff])
        value=coeff[-1];slope=finite(slope-h*(r0+r1)/2,'signed slope update')
    if not value.contains(0):raise ArithmeticError('Auxiliary endpoint misses zero')
    return lo,hi


def remainder_upper(nodes,delta):
    slope=arb(0)
    for a,b,d in zip(nodes[:-1],nodes[1:],delta):slope=finite(slope+d*(b-a)*(1-(a+b)/2),'remainder initial slope')
    value=arb(0);hi=arb(0)
    for a,b,d in zip(nodes[:-1],nodes[1:],delta):
        h=b-a;coeff=[value,value+slope*h/2,value+slope*h-d*h*h/2]
        for v in coeff:finite(v,'remainder Bernstein coefficient')
        hi=max(hi,*[v.abs_upper() for v in coeff])
        value=coeff[-1];slope=finite(slope-d*h,'remainder slope update')
    if not value.contains(0):raise ArithmeticError('Remainder endpoint misses zero')
    return hi


def certify(problem,n):
    if problem.boundary_error!=(0,0):raise ValueError('Need zero boundary error')
    start=time.perf_counter();problem.guard();nodes=make_mesh(n)
    r=[finite(problem.residual(x),'point residual') for x in nodes];delta=[]
    for a,b in zip(nodes[:-1],nodes[1:]):
        cell=a.lower().union(b.upper());K=finite(problem.r_second(cell),'residual second derivative').abs_upper();h=b-a
        delta.append(finite(K*h*h/8,'local remainder').upper())
    zl,zu=polynomial_norm(nodes,r);du=remainder_upper(nodes,delta)
    result=dict(lower=max(arb(0),(zl-du).lower()),upper=(zu+du).upper(),
                z_lower=zl,z_upper=zu,remainder_upper=du,intervals=n,bits=problem.bits,
                elapsed_s=time.perf_counter()-start)
    if exact_fraction(result['lower'])>exact_fraction(result['upper']):raise ArithmeticError('Reversed bounds')
    bounds(result)
    return result
