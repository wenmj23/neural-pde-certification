"""Exploratory outward interval backend for 1D piecewise-linear residual bounds.

Uses only basic interval arithmetic, exp and sin in mpmath.iv. The backend is
documented as experimental. This is not an independently audited certificate.
Binary-float network weights are treated as exact real constants. Mesh points
are the rational i/N, enclosed by interval division. Bernstein hulls bound the
continuous auxiliary polynomials without uncertain root finding.
"""
from pathlib import Path
import json,time
from mpmath import iv

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/generated/09_interval_prototype'


def zero():return iv.mpf(0)
def upper_abs(x):return abs(x).b


class SineResidual:
    def __init__(self,k):self.k=k;self.A=iv.mpf('0.0001');self.omega=k*iv.pi
    def residual(self,x):return self.A*self.omega**2*iv.sin(self.omega*x)
    def r_second(self,x):return -self.A*self.omega**4*iv.sin(self.omega*x)


class NeuralResidual:
    def __init__(self,weights,k=1):
        self.k=k;self.omega=k*iv.pi
        self.w=[iv.mpf(float(v[0])) for v in weights['hidden.weight']]
        self.b=[iv.mpf(float(v)) for v in weights['hidden.bias']]
        self.a=[iv.mpf(float(v)) for v in weights['output.weight'][0]]
        self.c=iv.mpf(float(weights['output.bias'][0]))
    def derivatives(self,x,order):
        vals=[self.c]+[zero() for _ in range(order)]
        for w,b,a in zip(self.w,self.b,self.a):
            t=1-2/(iv.exp(2*(w*x+b))+1)
            t2=t**2;t3=t**3
            ds=[t,1-t2,-2*t+2*t3]
            if order>=3:ds.append(-2+8*t2-6*t**4)
            if order>=4:ds.append(16*t-40*t3+24*t**5)
            for j in range(order+1):vals[j]+=a*w**j*ds[j]
        return vals
    def residual(self,x):
        n,n1,n2=self.derivatives(x,2)
        return -x*(1-x)*n2-2*(1-2*x)*n1+2*n-self.omega**2*iv.sin(self.omega*x)
    def r_second(self,x):
        _,_,n2,n3,n4=self.derivatives(x,4)
        return -x*(1-x)*n4-4*(1-2*x)*n3+12*n2+self.omega**4*iv.sin(self.omega*x)


def polynomial_enclosure(nodes,r):
    # Integrate the piecewise linear source; track all arithmetic as intervals.
    slope=zero()
    for a,b,r0,r1 in zip(nodes[:-1],nodes[1:],r[:-1],r[1:]):
        h=b-a
        slope+=h*(1-a)*(r0+r1)/2-h**2*(r0+2*r1)/6
    value=zero();lo=zero();hi=zero()
    for a,b,r0,r1 in zip(nodes[:-1],nodes[1:],r[:-1],r[1:]):
        h=b-a
        b0=value
        b1=value+slope*h/3
        b2=value+2*slope*h/3-r0*h**2/6
        b3=value+slope*h-h**2*(2*r0+r1)/6
        lo=max(lo,abs(b0).a,abs(b3).a)
        hi=max(hi,*[upper_abs(v) for v in (b0,b1,b2,b3)])
        value=b3;slope-=h*(r0+r1)/2
    assert 0 in value
    return lo,hi


def positive_remainder_upper(nodes,delta):
    slope=zero()
    for a,b,d in zip(nodes[:-1],nodes[1:],delta):
        h=b-a;slope+=d*h*(1-(a+b)/2)
    value=zero();hi=zero()
    for a,b,d in zip(nodes[:-1],nodes[1:],delta):
        h=b-a
        b0=value;b1=value+slope*h/2;b2=value+slope*h-d*h**2/2
        hi=max(hi,upper_abs(b0),upper_abs(b1),upper_abs(b2))
        value=b2;slope-=d*h
    assert 0 in value
    return hi


def run(problem,intervals):
    start=time.perf_counter()
    nodes=[iv.mpf(i)/intervals for i in range(intervals+1)]
    r=[problem.residual(x) for x in nodes]
    delta=[]
    for a,b in zip(nodes[:-1],nodes[1:]):
        cell=iv.mpf([a.a,b.b])
        K=upper_abs(problem.r_second(cell))
        delta.append((K*(b-a)**2/8).b)
    zl,zu=polynomial_enclosure(nodes,r)
    du=positive_remainder_upper(nodes,delta)
    lower=max(zero(),(zl-du).a);upper=(zu+du).b
    return dict(intervals=intervals,lower_interval=str(lower),upper_interval=str(upper),
                lower_display=float(lower),upper_display=float(upper),
                remainder_upper_interval=str(du),elapsed_s=time.perf_counter()-start),lower,upper


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    fixture=json.loads((ROOT/'configs/frozen_networks_round07.json').read_text())
    weights=next(row['weights'] for row in fixture['networks'] if row['label']=='PINN_k1_seed7')
    records=[]
    for precision,kind,k,n in [(30,'sine',1,128),(30,'sine',20,1000),(50,'sine',20,1000),
                                (30,'network',1,128),(30,'network',1,512),(50,'network',1,512)]:
        iv.dps=precision
        problem=SineResidual(k) if kind=='sine' else NeuralResidual(weights,k)
        row,lower,upper=run(problem,n)
        row.update(precision_digits=precision,kind=kind,k=k)
        if kind=='sine':
            reference=iv.mpf('0.0001');assert lower<=reference.a and upper>=reference.b
            row['upper_ratio_display']=float(upper)/1e-4
        else:
            reference=iv.mpf(0.00010329992292423817)
            assert lower<=reference and upper>=reference
            row['upper_ratio_to_reference_grid_display']=float(upper)/float(reference)
        records.append(row)
        (OUT/'records.json').write_text(json.dumps(records,indent=2))
        print(json.dumps(row),flush=True)
    (OUT/'scope.txt').write_text('Exploratory interval prototype using mpmath.iv 1.3.0, documented as experimental. Exact binary-float weights, rational mesh, interval residual and derivative evaluation, Bernstein polynomial hull, outward endpoints. No independent software audit; do not label this a production-grade universal verifier. Display floats are approximate; interval strings retain the reported enclosure.\n')


if __name__=='__main__':main()
