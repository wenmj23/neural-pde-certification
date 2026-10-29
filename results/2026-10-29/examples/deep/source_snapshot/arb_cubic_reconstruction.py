"""Experimental cubic residual interpolation and a moment remainder majorant.

Nodes 0,1/3,2/3,1 give max|nodal product|=1/81; delta=K4*h^4/1944.
Requires a residual fourth derivative bound (e.g. candidate C6, forcing C4).
"""
from fractions import Fraction as F
from math import comb
from flint import arb
from arb_deep_network import DeepPoisson,activate
from arb_reconstruction import Verifier,norm_enclosure,point
from arb_shaped_remainder import power_to_bernstein,split_exact
from arb_certification import finite
from arb_problem_families import SineError,BumpError


def multiply(a,b):
    result=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):result[i+j]+=x*y
    return result


NODES=[F(0),F(1,3),F(2,3),F(1)]
LAGRANGE=[]
for i,node in enumerate(NODES):
    p=[F(1)];denominator=F(1)
    for j,other in enumerate(NODES):
        if i!=j:p=multiply(p,[-other,F(1)]);denominator*=node-other
    LAGRANGE.append([v/denominator for v in p])
PRODUCT=[F(1)]
for node in NODES:PRODUCT=multiply(PRODUCT,[-node,F(1)])


def shift(poly,left,width):
    return [width**j*sum(poly[k]*comb(k,j)*left**(k-j) for k in range(j,len(poly))) for j in range(len(poly))]


def cubic_shape_constants():
    pieces=[];mass=F(0)
    for left,right,sign in [(F(0),F(1,3),-1),(F(1,3),F(2,3),1),(F(2,3),F(1),-1)]:
        source=shift([81*sign*x for x in PRODUCT],left,right-left)
        mass+=(right-left)*sum(c/F(j+1) for j,c in enumerate(source));pieces.append((right-left,source))
    value=slope=F(0);coefficient_upper=F(0)
    for width,source in pieces:
        q=[-v for v in source];q[0]+=mass
        power=[value,slope*width]+[width*width*c/F((j+1)*(j+2)) for j,c in enumerate(q)]
        subpieces=[power_to_bernstein(power)]
        for _ in range(4):subpieces=[child for p in subpieces for child in split_exact(p)]
        coefficient_upper=max(coefficient_upper,max(max(p) for p in subpieces))
        value=sum(power);slope+=width*sum(c/F(j+1) for j,c in enumerate(q))
    if value!=0 or slope!=0:raise ArithmeticError('Cubic moment identity endpoint mismatch')
    return mass,coefficient_upper


CUBIC_AVERAGE,CUBIC_CORRECTION=cubic_shape_constants()


class SixthOrderDeep(DeepPoisson):
    def __init__(self,fixture,label,sha,bits=103):super().__init__(fixture,label,sha,'chain',bits)
    def network_jet(self,x,order):
        if order<=5:return super().network_jet(x,order)
        if order!=6:raise ValueError('Cubic prototype supports at most six derivatives')
        self.guard();finite(x,'sixth-order input');self.forward_calls+=1
        state=[[x,arb(1)]+[arb(0)]*5]
        for layer in range(len(self.matrices)):
            state=self._affine(layer,state,6)
            if layer<len(self.matrices)-1:state=[activate(z,'chain') for z in state]
        return state[0]
    def r_fourth(self,x):
        n=self.network_jet(x,6)
        return finite(-720*(self.q(x)*n[6]+(1-2*x)*n[5]-n[4])-self.wp[6]*(self.omega*x).sin(),'fourth residual derivative')


class SineSix(SineError):
    def r_fourth(self,x):
        self.guard();result=arb(0)
        for w,c in zip(self.omegas,self.co):result+=self.A*c*w*w*w*w*w*w*(w*x).sin()
        return finite(result,'analytic fourth residual derivative')


class BumpSix(BumpError):
    def __init__(self,width,bits=103):
        super().__init__(width,bits);self.coefficients.append([-120,0,720,0,-480,0,64])
    def r_fourth(self,x):return finite(-self._e_derivative(x,6),'bump fourth residual derivative')


class CubicVerifier(Verifier):
    def __init__(self,problem,mode='moment',subdivision=2,cache=True):
        if mode not in ('constant','moment'):raise ValueError('Unknown cubic remainder mode')
        if not hasattr(problem,'r_fourth'):raise ValueError('Fourth residual derivative bound required')
        super().__init__(problem,degree=2,subdivision=subdivision,cache=cache)
        self.degree=3;self.mode=mode
        self.average=arb(CUBIC_AVERAGE.numerator)/CUBIC_AVERAGE.denominator
        self.correction=arb(CUBIC_CORRECTION.numerator)/CUBIC_CORRECTION.denominator
    def prepare(self,nodes):
        self.guard()
        if len(nodes)<2 or nodes[0]!=0 or nodes[-1]!=1 or any(a>=b for a,b in zip(nodes[:-1],nodes[1:])):raise ValueError('Invalid rational mesh')
        positions=set(nodes)
        for a,b in zip(nodes[:-1],nodes[1:]):positions.update(a+(b-a)*t for t in NODES)
        values={x:self._residual(x) for x in sorted(positions)};cells=[]
        for a,b in zip(nodes[:-1],nodes[1:]):
            if self.cache and (a,b) in self.cell_cache:cells.append(self.cell_cache[(a,b)]);continue
            aa,bb=point(a),point(b);h=finite(bb-aa,'cubic cell width')
            if h.lower()<=0:raise ValueError('Unseparated mesh')
            source=[values[a+(b-a)*t] for t in NODES]
            power=[]
            for j in range(4):
                p=arb(0)
                for i in range(4):p+=source[i]*LAGRANGE[i][j].numerator/LAGRANGE[i][j].denominator
                power.append(finite(p,'cubic residual coefficient'))
            K=finite(self.problem.r_fourth(aa.lower().union(bb.upper())),'fourth residual derivative').abs_upper()
            delta=finite(K*h*h*h*h/1944,'cubic remainder').upper();self.curvature_evaluations+=1
            item=dict(a=aa,b=bb,h=h,p=power,delta=delta)
            if self.cache:self.cell_cache[(a,b)]=item
            cells.append(item)
        return cells
    def _norms(self,cells):
        self.polynomial_cell_visits+=len(cells);slope=arb(0);dslope=arb(0)
        for c in cells:
            h,p=c['h'],c['p'];mass=sum((v/(j+1) for j,v in enumerate(p)),arb(0));first=sum((v/(j+2) for j,v in enumerate(p)),arb(0))
            slope=finite(slope+h*((1-c['a'])*mass-h*first),'cubic initial slope')
            dslope=finite(dslope+c['delta']*h*(1-(c['a']+c['b'])/2),'cubic remainder initial slope')
        value=dvalue=arb(0);zl=zu=du=arb(0)
        for c in cells:
            h,p,delta=c['h'],c['p'],c['delta']
            power=[value,slope*h]+[-h*h*v/((j+1)*(j+2)) for j,v in enumerate(p)]
            coeff=[finite(v,'cubic auxiliary Bernstein coefficient') for v in power_to_bernstein(power)]
            low,high=norm_enclosure(coeff,self.depth);zl=max(zl,low);zu=max(zu,high)
            value=coeff[-1];slope=finite(slope-h*sum((v/(j+1) for j,v in enumerate(p)),arb(0)),'cubic slope')
            dc=[dvalue,dvalue+dslope*h/2,dvalue+dslope*h-delta*h*h/2]
            _,high=norm_enclosure(dc,self.depth);du=max(du,high)
            dvalue=dc[-1];dslope=finite(dslope-delta*h,'cubic remainder slope')
        if not value.contains(0) or not dvalue.contains(0):raise ArithmeticError('Cubic integrated endpoint misses zero')
        if self.mode=='moment':
            correction=max(finite(c['delta']*c['h']*c['h'],'cubic moment correction').upper() for c in cells)
            candidate=finite(self.average*du+self.correction*correction,'cubic moment bound').upper()
            du=min(du,candidate)
        return zl,zu,du
