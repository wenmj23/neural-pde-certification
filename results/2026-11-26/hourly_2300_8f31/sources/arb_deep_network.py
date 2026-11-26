"""Validated derivative jets for frozen dense tanh networks on one input.

Coefficient j encloses N^(j)(x)/j! for every x in an input interval.
Methods: exp (native formal exp/reciprocal), chain (factored tanh derivatives
and native formal composition), centered (chain plus mean-value contraction).
Jets are derivative enclosures, not a truncated approximation without remainder.
"""
from fractions import Fraction
from math import factorial,isfinite
import hashlib
import json
from pathlib import Path
from flint import arb,arb_mat,arb_series,ctx
from arb_certification import finite,exact_fraction,square,int_power
from arb_problem_families import SmoothProblem,derivative_polys


def activation_polys(t,order):
    result=derivative_polys(t,min(order,5))
    if order>=6:
        t2=square(t)
        result.append(-16*t*(1-t2)*(17+t2*(-60+45*t2)))
    return [finite(x,'tanh derivative polynomial') for x in result]


def activate(z,method='chain'):
    order=len(z)-1
    if method=='exp':
        s=arb_series(z,prec=order+1)
        y=1-2/((2*s).exp()+1)
        return [finite(y[j],'native exponential jet') for j in range(order+1)]
    t=z[0].tanh()
    if not t.is_finite():
        # On a tested wide finite interval the binding returned indeterminate.
        # tanh is increasing, so endpoint bounds are an independent safe form.
        low=finite(z[0].lower().tanh(),'lower endpoint tanh').lower()
        high=finite(z[0].upper().tanh(),'upper endpoint tanh').upper()
        t=low.union(high)
    t=finite(t,'tanh value')
    if order==0:return [t]
    polys=activation_polys(t,order)
    delta=arb_series([arb(0)]+z[1:],prec=order+1)
    power=arb_series([1],prec=order+1);result=arb_series([t],prec=order+1)
    for j in range(1,order+1):
        power=power*delta
        result=result+power*(polys[j]/factorial(j))
    return [finite(result[j],'factored composition jet') for j in range(order+1)]


def contract(jet,center,radius):
    values=list(jet)
    for j in range(len(jet)-2,-1,-1):
        error=finite(radius*(j+1)*values[j+1].abs_upper(),'mean value radius').upper()
        candidate=finite(center[j]+arb(0,error),'mean value enclosure')
        values[j]=finite(values[j].intersection(candidate),'mean value intersection')
    return values


class DeepPoisson(SmoothProblem):
    def __init__(self,fixture,label,expected_sha256,method='chain',bits=103):
        self.setup(bits)
        if ctx.cap<7:raise ValueError('Series truncation cap must be at least seven')
        self.cap=ctx.cap
        if method not in ('exp','chain','centered'):raise ValueError('Unknown derivative method')
        data=Path(fixture).read_bytes()
        if hashlib.sha256(data).hexdigest()!=expected_sha256:raise ValueError('Frozen depth fixture changed')
        matches=[r for r in json.loads(data)['networks'] if r['label']==label]
        if len(matches)!=1:raise ValueError('Unknown/ambiguous deep network')
        self.raw=matches[0];self.name=label;self.method=method;self.frequency=self.raw['frequency']
        if type(self.frequency) is not int or self.frequency<1:raise ValueError('Invalid frequency')
        self.omega=self.frequency*arb.pi();self.wp=[int_power(self.omega,j) for j in range(7)]
        def value(v):
            if isinstance(v,bool) or not isinstance(v,(int,float)) or not isfinite(v):raise ValueError('Invalid frozen weight')
            n,d=float(v).as_integer_ratio();a=arb(n)/d
            if exact_fraction(a)!=Fraction(n,d):raise ArithmeticError('Frozen binary64 value not imported exactly')
            return a
        self.matrices=[];self.biases=[];previous=1
        layers=self.raw['layers']
        if not 2<=len(layers)<=8:raise ValueError('Unsupported network depth')
        for i,layer in enumerate(layers):
            weights=layer['weight'];bias=layer['bias'];expected='identity' if i==len(layers)-1 else 'tanh'
            if layer['activation']!=expected:raise ValueError('Unsupported activation')
            if not weights or len(weights)!=len(bias) or any(len(row)!=previous for row in weights):raise ValueError('Wrong affine layer shape')
            if i==len(layers)-1 and len(weights)!=1:raise ValueError('Scalar output required')
            self.matrices.append(arb_mat([[value(x) for x in row] for row in weights]))
            self.biases.append([value(x) for x in bias]);previous=len(weights)
        self.forward_calls=0

    def guard(self):
        super().guard()
        if ctx.cap!=self.cap:raise ValueError('Series truncation cap changed; rebuild problem/cache')

    def _affine(self,layer,state,order):
        result=self.matrices[layer]*arb_mat(state)
        return [[finite(result[i,j]+(self.biases[layer][i] if j==0 else 0),'affine derivative jet')
                 for j in range(order+1)] for i in range(len(self.biases[layer]))]

    def network_jet(self,x,order):
        self.guard();finite(x,'network input')
        if type(order) is not int or not 0<=order<=5:raise ValueError('Supported network derivative orders are zero through five')
        self.forward_calls+=1
        centered=self.method=='centered' and x.rad().upper()>arb(1)/10**20
        work_order=order+1 if centered else order
        state=[[x]+([arb(1)]+[arb(0)]*(work_order-1) if work_order else [])]
        if centered:
            centers=[[x.mid()]+([arb(1)]+[arb(0)]*(work_order-1) if work_order else [])]
            radius=x.rad().upper()
        for layer in range(len(self.matrices)):
            affine=self._affine(layer,state,work_order)
            if centered:
                ca=self._affine(layer,centers,work_order)
                affine=[contract(z,c,radius) for z,c in zip(affine,ca)]
            if layer<len(self.matrices)-1:
                method='exp' if self.method=='exp' else 'chain'
                state=[activate(z,method) for z in affine]
                if centered:
                    centers=[activate(z,'chain') for z in ca]
                    state=[contract(z,c,radius) for z,c in zip(state,centers)]
            else:state=affine
        return state[0][:order+1]

    def derivatives(self,x,orders):
        jet=self.network_jet(x,max(orders))
        return {j:finite(jet[j]*factorial(j),'network derivative') for j in orders}

    def solution(self,x):return finite(self.q(x)*self.network_jet(x,0)[0],'network solution')
    def error(self,x):return self.solution(x)-(self.omega*x).sin()
    def residual(self,x):
        n=self.network_jet(x,2)
        return finite(-2*(self.q(x)*n[2]+(1-2*x)*n[1]-n[0])-self.wp[2]*(self.omega*x).sin(),'residual')
    def r_second(self,x):
        n=self.network_jet(x,4)
        return finite(-24*(self.q(x)*n[4]+(1-2*x)*n[3]-n[2])+self.wp[4]*(self.omega*x).sin(),'second residual derivative')
    def r_third(self,x):
        n=self.network_jet(x,5)
        return finite(-120*(self.q(x)*n[5]+(1-2*x)*n[4]-n[3])+self.wp[5]*(self.omega*x).cos(),'third residual derivative')
