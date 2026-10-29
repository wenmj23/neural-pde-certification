"""Position-dependent Lagrange remainder for quadratic residual reconstruction.

For a cell t=(x-a)/h, the old constant bound is delta. A safe shaped
source is 12*sqrt(3)*delta*|t(t-1/2)(t-1)|. Integrate it as two cubic
pieces, or use its exact moments and a polynomial correction majorant.
The frozen-data/precision and finite-arithmetic guards of Verifier remain.
"""
from fractions import Fraction
from math import comb, ceil, isfinite
from flint import arb
from arb_certification import finite, exact_fraction
from arb_reconstruction import Verifier, norm_enclosure


def power_to_bernstein(power):
    degree = len(power)-1
    return [sum((power[j]*comb(k,j)/comb(degree,j) for j in range(k+1)), power[0]*0)
            for k in range(degree+1)]


def split_exact(coeff):
    row=list(coeff);left=[row[0]];right=[row[-1]]
    while len(row)>1:
        row=[(a+b)/2 for a,b in zip(row[:-1],row[1:])]
        left.append(row[0]);right.append(row[-1])
    return left,list(reversed(right))


def correction_constant(depth=4):
    # J(t)/sqrt(3) on t=u/2. Right half is its reflection.
    power=[Fraction(0),Fraction(0),Fraction(3,64),Fraction(-1,8),Fraction(3,32),Fraction(-3,160)]
    pieces=[power_to_bernstein(power)]
    for _ in range(depth):
        pieces=[child for p in pieces for child in split_exact(p)]
    return max(max(p) for p in pieces)


KAPPA = correction_constant()


def arb_fraction(q):
    return arb(q.numerator)/q.denominator


def shaped_norm(cells, depth=2):
    """Continuous norm bracket for -D''=shaped_source, zero Dirichlet data."""
    sqrt3=arb(3).sqrt();average=3*sqrt3/8
    slope=arb(0)
    for c in cells:
        slope=finite(slope+average*c['delta']*c['h']*(1-(c['a']+c['b'])/2), 'shaped initial slope')
    value=arb(0);lower=arb(0);upper=arb(0)
    for c in cells:
        amp=12*sqrt3*c['delta'];width=c['h']/2
        halves=[[arb(0),amp/4,-3*amp/8,amp/8],
                [arb(0),amp/8,arb(0),-amp/8]]
        for source in halves:
            power=[value,slope*width]+[-width*width*source[j]/((j+1)*(j+2)) for j in range(4)]
            coeff=[finite(x,'quintic remainder coefficient') for x in power_to_bernstein(power)]
            lo,hi=norm_enclosure(coeff,depth);lower=max(lower,lo);upper=max(upper,hi)
            value=coeff[-1]
            mass=width*sum((source[j]/(j+1) for j in range(4)),arb(0))
            slope=finite(slope-mass,'shaped slope update')
    if not value.contains(0):
        raise ArithmeticError('Shaped endpoint misses zero')
    return lower,finite(upper,'shaped continuous norm')


class ShapedVerifier(Verifier):
    def __init__(self,problem,mode='moment',indicator='global',subdivision=2,cache=True):
        if mode not in ('constant','full','moment'):
            raise ValueError('Unknown remainder mode')
        if indicator not in ('global','peak'):
            raise ValueError('Unknown refinement indicator')
        super().__init__(problem,degree=2,subdivision=subdivision,cache=cache)
        self.mode=mode;self.indicator=indicator
        self.average=3*arb(3).sqrt()/8
        self.correction=arb_fraction(KAPPA)*arb(3).sqrt()
        self.last_peak=Fraction(1,2)

    def _norms(self,cells):
        zl,zu,constant=super()._norms(cells)
        self.constant_remainder=constant
        if self.mode=='constant':
            candidate=constant
        elif self.mode=='full':
            _,candidate=shaped_norm(cells,self.depth)
        else:
            correction=max(finite(c['delta']*c['h']*c['h'],'moment correction').upper() for c in cells)
            candidate=finite(self.average*constant+self.correction*correction,'moment remainder norm').upper()
        self.raw_remainder=finite(candidate,'raw shaped remainder norm').upper()
        # Both are independently valid bounds; retaining the old one prevents
        # numerical hull overestimation from worsening the final certificate.
        return zl,zu,min(constant,self.raw_remainder)

    def certify(self,nodes=None):
        r=super().certify(nodes)
        r.update(mode=self.mode,indicator=self.indicator,
                 constant_remainder_upper=self.constant_remainder,
                 raw_remainder_upper=self.raw_remainder,
                 indicator_peak=str(self.last_peak))
        return r

    def _peak(self,cells):
        # Mesh heuristic only: approximate the peak of the old concave D.
        # Floating arithmetic here never enters the final certificate.
        heights=[float(c['delta']) for c in cells]
        if not all(isfinite(x) and x>=0 for x in heights):
            raise ArithmeticError('Nonfinite mesh heuristic')
        slope=sum(d*float(c['h'])*(1-(float(c['a'])+float(c['b']))/2) for d,c in zip(heights,cells))
        if not isfinite(slope):
            raise ArithmeticError('Overflow in mesh heuristic')
        for d,c in zip(heights,cells):
            a,h=float(c['a']),float(c['h'])
            if d>0 and slope-d*h<=0:
                x=min(1.,max(0.,a+slope/d))
                return Fraction(x).limit_denominator(2**24)
            slope-=d*h
        return Fraction(1,2)

    def refine_to(self,n):
        if self.indicator=='global':
            return super().refine_to(n)
        if type(n) is not int or n<len(self.nodes)-1:
            raise ValueError('Invalid adaptive budget')
        self.guard()
        average=Fraction(1) if self.mode=='constant' else exact_fraction(self.average.upper())
        while len(self.nodes)-1<n:
            cells=self.prepare(self.nodes);peak=self._peak(cells);self.last_peak=peak;priorities=[]
            for i,c in enumerate(cells):
                a,b=self.nodes[i],self.nodes[i+1];h=b-a;mid=(a+b)/2;delta=exact_fraction(c['delta'])
                if b<=peak:
                    weight=average*h*(1-peak)*mid
                elif a>=peak:
                    weight=average*h*peak*(1-mid)
                else:
                    constant_weight=((1-peak)*(peak*peak-a*a)+peak*((1-peak)**2-(1-b)**2))/2
                    weight=constant_weight if self.mode=='constant' else min(constant_weight,average*h*peak*(1-peak))
                if weight<0:
                    raise ArithmeticError('Negative Green indicator')
                priorities.append((delta*weight,i))
            count=min(n-(len(self.nodes)-1),max(1,ceil(.25*(len(self.nodes)-1))))
            chosen=[i for _,i in sorted(priorities)[-count:]]
            self.nodes=sorted(self.nodes+[(self.nodes[i]+self.nodes[i+1])/2 for i in chosen])
