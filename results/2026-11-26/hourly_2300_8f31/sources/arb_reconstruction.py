"""Linear/quadratic signed residual reconstruction on arbitrary rational meshes.

Quadratic interpolation uses endpoints and midpoint with remainder
K3*h^3/(72*sqrt(3)). Its auxiliary solution is piecewise quartic. Bernstein
subdivision bounds continuous polynomial norms without new network evaluations.
All results remain conditional on the smooth model and interval backend.
"""
from fractions import Fraction
import math
import time
from flint import arb,ctx
from arb_certification import finite,exact_fraction,bounds


def point(x):return arb(x.numerator)/x.denominator


def split_bernstein(coeff):
    row=list(coeff);left=[row[0]];right=[row[-1]]
    while len(row)>1:
        row=[finite((a+b)/2,'Bernstein subdivision') for a,b in zip(row[:-1],row[1:])]
        left.append(row[0]);right.append(row[-1])
    return left,list(reversed(right))


def norm_enclosure(coeff,depth=0):
    if type(depth) is not int or not 0<=depth<=4:raise ValueError('Unsupported subdivision depth')
    for v in coeff:finite(v,'Bernstein coefficient')
    pieces=[list(coeff)]
    for _ in range(depth):pieces=[child for p in pieces for child in split_bernstein(p)]
    lo=arb(0);hi=arb(0)
    for p in pieces:
        lo=max(lo,p[0].abs_lower(),p[-1].abs_lower())
        hi=max(hi,*[c.abs_upper() for c in p])
    return lo,hi


class Verifier:
    def __init__(self,problem,degree=1,subdivision=0,cache=True):
        if degree not in (1,2):raise ValueError('Residual degree must be one or two')
        if problem.boundary_error!=(0,0):raise ValueError('Exactly zero error boundary required')
        if type(subdivision) is not int or not 0<=subdivision<=4:raise ValueError('Invalid subdivision depth')
        self.problem=problem;self.degree=degree;self.depth=subdivision;self.cache=cache
        self.node_cache={};self.cell_cache={};self.residual_evaluations=0;self.curvature_evaluations=0
        self.polynomial_cell_visits=0;self.nodes=[Fraction(i,16) for i in range(17)]
        problem.guard();self.bits=ctx.prec;self.quadratic_denominator=72*arb(3).sqrt()
    def guard(self):
        self.problem.guard()
        if ctx.prec!=self.bits:raise ValueError('Precision changed; rebuild cache')
    def _residual(self,x):
        if self.cache and x in self.node_cache:return self.node_cache[x]
        value=finite(self.problem.residual(point(x)),'cached residual');self.residual_evaluations+=1
        if self.cache:self.node_cache[x]=value
        return value
    def prepare(self,nodes):
        self.guard()
        if len(nodes)<2 or nodes[0]!=0 or nodes[-1]!=1 or any(a>=b for a,b in zip(nodes[:-1],nodes[1:])):
            raise ValueError('Strictly increasing rational mesh from zero to one required')
        positions=set(nodes)
        if self.degree==2:positions.update((a+b)/2 for a,b in zip(nodes[:-1],nodes[1:]))
        local={x:self._residual(x) for x in sorted(positions)}
        cells=[]
        for a,b in zip(nodes[:-1],nodes[1:]):
            key=(a,b)
            if self.cache and key in self.cell_cache:
                cells.append(self.cell_cache[key]);continue
            aa,bb=point(a),point(b);h=finite(bb-aa,'cell width')
            if h.lower()<=0:raise ValueError('Unseparated mesh')
            cell=aa.lower().union(bb.upper())
            derivative=self.problem.r_second(cell) if self.degree==1 else self.problem.r_third(cell)
            K=finite(derivative,'residual derivative bound').abs_upper();self.curvature_evaluations+=1
            delta=K*h*h/8 if self.degree==1 else K*h*h*h/self.quadratic_denominator
            item=dict(a=aa,b=bb,h=h,r0=local[a],r1=local[b],delta=finite(delta,'interpolation remainder').upper())
            if self.degree==2:item['rm']=local[(a+b)/2]
            if self.cache:self.cell_cache[key]=item
            cells.append(item)
        return cells
    def _norms(self,cells):
        self.polynomial_cell_visits+=len(cells);slope=arb(0);dslope=arb(0)
        for c in cells:
            a,b,h,r0,r1=c['a'],c['b'],c['h'],c['r0'],c['r1']
            if self.degree==1:term=h*(1-a)*(r0+r1)/2-h*h*(r0+2*r1)/6
            else:
                mid=(a+b)/2
                term=h*((1-a)*r0+4*(1-mid)*c['rm']+(1-b)*r1)/6
            slope=finite(slope+term,'auxiliary initial slope')
            dslope=finite(dslope+c['delta']*h*(1-(a+b)/2),'remainder initial slope')
        value=arb(0);dvalue=arb(0);zl=arb(0);zu=arb(0);du=arb(0)
        for c in cells:
            h,r0,r1,delta=c['h'],c['r0'],c['r1'],c['delta'];h2=h*h
            if self.degree==1:
                coeff=[value,value+slope*h/3,value+2*slope*h/3-r0*h2/6,
                       value+slope*h-h2*(2*r0+r1)/6]
                mass=h*(r0+r1)/2
            else:
                B=4*c['rm']-3*r0-r1;C=2*(r0+r1-2*c['rm'])
                coeff=[value,value+slope*h/4,value+slope*h/2-r0*h2/12,
                       value+3*slope*h/4-r0*h2/4-B*h2/24,
                       value+slope*h-r0*h2/2-B*h2/6-C*h2/12]
                mass=h*(r0+4*c['rm']+r1)/6
            low,high=norm_enclosure(coeff,self.depth);zl=max(zl,low);zu=max(zu,high)
            value=coeff[-1];slope=finite(slope-mass,'auxiliary slope update')
            dc=[dvalue,dvalue+dslope*h/2,dvalue+dslope*h-delta*h2/2]
            _,high=norm_enclosure(dc,self.depth);du=max(du,high)
            dvalue=dc[-1];dslope=finite(dslope-delta*h,'remainder slope update')
        if not value.contains(0) or not dvalue.contains(0):raise ArithmeticError('Integrated endpoint misses zero')
        return zl,zu,du
    def certify(self,nodes=None):
        start=time.perf_counter();self.guard()
        nodes=self.nodes if nodes is None else list(nodes)
        cells=self.prepare(nodes);zl,zu,du=self._norms(cells)
        result=dict(lower=max(arb(0),(zl-du).lower()),upper=(zu+du).upper(),
            z_lower=zl,z_upper=zu,remainder_upper=du,intervals=len(nodes)-1,
            degree=self.degree,subdivision=self.depth,elapsed_s=time.perf_counter()-start,
            residual_evaluations=self.residual_evaluations,curvature_evaluations=self.curvature_evaluations,
            polynomial_cell_visits=self.polynomial_cell_visits)
        if exact_fraction(result['lower'])>exact_fraction(result['upper']):raise ArithmeticError('Reversed bounds')
        bounds(result);return result
    def uniform(self,n):
        if type(n) is not int or n<1:raise ValueError('Positive integer mesh budget required')
        self.nodes=[Fraction(i,n) for i in range(n+1)]
    def refine_to(self,n):
        if type(n) is not int or n<len(self.nodes)-1:raise ValueError('Invalid adaptive budget')
        self.guard()
        while len(self.nodes)-1<n:
            cells=self.prepare(self.nodes);priorities=[]
            for i,c in enumerate(cells):
                a,b=self.nodes[i],self.nodes[i+1]
                if a<=Fraction(1,2)<=b:weight=arb(1)/4
                else:weight=max((point(a)*(1-point(a))).upper(),(point(b)*(1-point(b))).upper())
                score=finite(c['delta']*c['h']*weight,'adaptive indicator').upper()
                priorities.append((exact_fraction(score),i))
            count=min(n-(len(self.nodes)-1),max(1,math.ceil(.25*(len(self.nodes)-1))))
            selected=[i for _,i in sorted(priorities)[-count:]]
            self.nodes=sorted(self.nodes+[(self.nodes[i]+self.nodes[i+1])/2 for i in selected])
