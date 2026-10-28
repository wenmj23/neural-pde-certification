"""Working 1D Poisson certificate algorithms (ordinary float64, not outward-rounded).

Independent of the solver: a problem supplies residual values and valid local
bounds for residual derivatives. The caller must ensure zero error at boundaries.
"""
from dataclasses import dataclass
from functools import lru_cache
import time
import numpy as np
from numpy.polynomial import polynomial as pol


def interval_product(a,b):
    lo,hi=a;ll,uu=b
    v=np.array([lo*ll,lo*uu,hi*ll,hi*uu])
    return v.min(axis=0),v.max(axis=0)


@lru_cache(maxsize=None)
def critical_points(coeff):
    dc=pol.polyder(coeff)
    return tuple(z.real for z in (pol.polyroots(dc) if len(dc)>1 else []) if abs(z.imag)<=1e-12)


def polynomial_range(lo,hi,coeff):
    a,b=pol.polyval(lo,coeff),pol.polyval(hi,coeff)
    lower,upper=np.minimum(a,b),np.maximum(a,b)
    for z in critical_points(tuple(coeff)):
        mask=(lo<=z)&(z<=hi);v=pol.polyval(z,coeff)
        lower=np.where(mask,np.minimum(lower,v),lower)
        upper=np.where(mask,np.maximum(upper,v),upper)
    return lower,upper


def trig_range(left,right,frequency,kind='cos'):
    lo,hi=frequency*left,frequency*right
    if kind=='sin':lo=lo-np.pi/2;hi=hi-np.pi/2
    a,b=np.cos(lo),np.cos(hi)
    lower,upper=np.minimum(a,b),np.maximum(a,b)
    maxima=np.ceil(lo/(2*np.pi))<=np.floor(hi/(2*np.pi))
    minima=np.ceil((lo-np.pi)/(2*np.pi))<=np.floor((hi-np.pi)/(2*np.pi))
    return np.where(minima,-1.,lower),np.where(maxima,1.,upper)


def q_range(left,right):
    lo=np.minimum(left*(1-left),right*(1-right))
    hi=np.maximum(left*(1-left),right*(1-right))
    return lo,np.where((left<=.5)&(right>=.5),.25,hi)


def integrate_piecewise_constant(nodes,forcing):
    h=np.diff(nodes);mid=(nodes[:-1]+nodes[1:])/2
    mass=h*forcing;moment=mass*mid
    prefix_mass=np.r_[0.,np.cumsum(mass)]
    prefix_moment=np.r_[0.,np.cumsum(moment)]
    slope0=np.sum(mass*(1-mid))
    values=slope0*nodes-nodes*prefix_mass+prefix_moment
    slopes=slope0-prefix_mass[:-1]
    potential=(forcing!=0)&(np.signbit(slopes)==np.signbit(forcing))&(abs(slopes)<=abs(forcing)*h)
    root=np.divide(slopes,forcing,out=np.full_like(h,np.inf),where=potential)
    valid=(root>0)&(root<h)
    t=np.where(valid,root,0.)
    extreme=values[:-1]+slopes*t-.5*forcing*t*t
    maximum=max(np.max(abs(values)),np.max(np.where(valid,abs(extreme),0.)))
    return float(maximum),values


def integrate_piecewise_linear(nodes,residual):
    h=np.diff(nodes);r0,r1=residual[:-1],residual[1:];dr=r1-r0
    mass=h*(r0+r1)/2
    moment=(nodes[:-1]+h/2)*mass+h*h*dr/12
    pm=np.r_[0.,np.cumsum(mass)];pq=np.r_[0.,np.cumsum(moment)]
    slope0=np.sum(mass*(1-nodes[:-1]-h/2)-h*h*dr/12)
    values=slope0*nodes-nodes*pm+pq;slopes=slope0-pm[:-1]
    # Roots of (dr/2h)*t^2+r0*t-slopes=0, using a stable quadratic formula.
    aa=dr/(2*h);bb=r0;disc=bb*bb+4*aa*slopes
    sq=np.sqrt(np.maximum(disc,0.))
    qq=-.5*(bb+np.where(bb>=0,1.,-1.)*sq)
    roots=[np.divide(qq,aa,out=np.full_like(h,np.inf),where=aa!=0),
           np.divide(-slopes,qq,out=np.full_like(h,np.inf),where=qq!=0),
           np.divide(slopes,bb,out=np.full_like(h,np.inf),where=(aa==0)&(bb!=0))]
    maximum=float(np.max(abs(values)))
    for root in roots:
        valid=(root>0)&(root<h)&(disc>=0)
        t=np.where(valid,root,0.)
        val=values[:-1]+slopes*t-r0*t*t/2-dr*t**3/(6*h)
        maximum=max(maximum,float(np.max(np.where(valid,abs(val),0.))))
    return maximum,values


def certificate(problem,nodes,method='linear'):
    nodes=np.asarray(nodes,dtype=float)
    if method not in ('constant','linear'):raise ValueError('Unknown reconstruction method')
    if nodes.ndim!=1 or len(nodes)<2 or nodes[0]!=0 or nodes[-1]!=1 or np.any(np.diff(nodes)<=0):
        raise ValueError('Nodes must increase strictly from 0 to 1')
    start=time.perf_counter();h=np.diff(nodes)
    r=problem.residual(nodes)
    order=1 if method=='constant' else 2
    bound=problem.derivative_bound(nodes,order)
    if method=='constant':
        delta=bound*h/2
        Z,z=integrate_piecewise_constant(nodes,(r[:-1]+r[1:])/2)
    else:
        delta=bound*h*h/8
        Z,z=integrate_piecewise_linear(nodes,r)
    D,d=integrate_piecewise_constant(nodes,delta)
    return dict(lower=max(0.,Z-D),upper=Z+D,Z=Z,D=D,
                elapsed_ms=(time.perf_counter()-start)*1000,
                nodes=len(nodes),z=z,d=d,derivative_bound=bound)


def adaptive_nodes(problem,intervals,method='linear',initial_nodes=None):
    """Batch refinement by an upper bound on each cell's Green-weighted remainder.

    Reports all interval bound evaluations, not only the final number of cells.
    No exact solution/error is used to select the mesh.
    """
    start=time.perf_counter()
    nodes=np.linspace(0,1,min(16,intervals)+1) if initial_nodes is None else np.array(initial_nodes,copy=True)
    if len(nodes)-1>intervals:raise ValueError('Target mesh smaller than current mesh')
    queries=0;order=1 if method=='constant' else 2
    while len(nodes)-1<intervals:
        h=np.diff(nodes);K=problem.derivative_bound(nodes,order);queries+=len(h)
        delta=K*h/2 if order==1 else K*h*h/8
        _,weight=q_range(nodes[:-1],nodes[1:])
        indicator=delta*h*weight
        count=min(intervals-len(h),max(1,int(np.ceil(.3*len(h)))))
        chosen=np.argsort(indicator,kind='stable')[-count:]
        nodes=np.sort(np.r_[nodes,(nodes[chosen]+nodes[chosen+1])/2])
    return nodes,dict(mesh_ms=(time.perf_counter()-start)*1000,bound_evaluations=queries)


@dataclass
class SineError:
    frequency:int
    amplitude:float=1e-4

    @property
    def name(self):return f'sine_{self.frequency}'
    @property
    def exact_error(self):return abs(self.amplitude)
    def error(self,x):return self.amplitude*np.sin(self.frequency*np.pi*x)
    def residual(self,x):return (self.frequency*np.pi)**2*self.error(x)
    def derivative_bound(self,nodes,order):
        kind='cos' if order==1 else 'sin'
        lo,hi=trig_range(nodes[:-1],nodes[1:],self.frequency*np.pi,kind)
        return abs(self.amplitude)*(self.frequency*np.pi)**(order+2)*np.maximum(abs(lo),abs(hi))


@dataclass
class BumpError:
    width:float
    amplitude:float=1e-4

    @property
    def name(self):return f'bump_{self.width}'
    @property
    def exact_error(self):return abs(self.amplitude)
    def error(self,x):return 4*self.amplitude*x*(1-x)*np.exp(-((x-.5)/self.width)**2)
    def residual(self,x):
        t=(x-.5)/self.width;g=np.exp(-t*t);q=x*(1-x);qp=1-2*x
        g1=-2*t/self.width*g;g2=(4*t*t-2)/self.width**2*g
        return -4*self.amplitude*(q*g2+2*qp*g1-2*g)
    def derivative_bound(self,nodes,order):
        left,right=nodes[:-1],nodes[1:];lo=(left-.5)/self.width;hi=(right-.5)/self.width
        tmin=np.where((lo<=0)&(hi>=0),0.,np.minimum(abs(lo),abs(hi)))
        tmax=np.maximum(abs(lo),abs(hi));gr=(np.exp(-tmax*tmax),np.exp(-tmin*tmin))
        polys=[[1.],[0.,-2.],[-2.,0.,4.],[0.,12.,0.,-8.],[12.,0.,-48.,0.,16.]]
        gs=[]
        for j in range(5):
            a,b=interval_product(polynomial_range(lo,hi,polys[j]),gr)
            gs.append((a/self.width**j,b/self.width**j))
        qr=q_range(left,right);qpr=(1-2*right,1-2*left)
        if order==1:
            p=interval_product(qr,gs[3]);s=interval_product(qpr,gs[2])
            lower=p[0]+3*s[0]-6*gs[1][1];upper=p[1]+3*s[1]-6*gs[1][0]
        else:
            p=interval_product(qr,gs[4]);s=interval_product(qpr,gs[3])
            lower=p[0]+4*s[0]-12*gs[2][1];upper=p[1]+4*s[1]-12*gs[2][0]
        return 4*abs(self.amplitude)*np.maximum(abs(lower),abs(upper))


@dataclass
class NeuralResidual:
    params:tuple
    frequency:int
    label:str='network'

    @property
    def name(self):return self.label
    def residual(self,x):
        w,b,a,c=self.params;t=np.tanh(np.asarray(x)[...,None]*w+b)
        n=c+(t*a).sum(axis=-1)
        n1=((1-t*t)*(a*w)).sum(axis=-1)
        n2=((-2*t+2*t**3)*(a*w**2)).sum(axis=-1)
        return -(x*(1-x)*n2+2*(1-2*x)*n1-2*n)-(self.frequency*np.pi)**2*np.sin(self.frequency*np.pi*x)
    def derivative_bound(self,nodes,order):
        w,b,a,_=self.params;left,right=nodes[:-1],nodes[1:]
        z1,z2=left[:,None]*w+b,right[:,None]*w+b
        lo,hi=np.tanh(np.minimum(z1,z2)),np.tanh(np.maximum(z1,z2))
        polys={1:[1.,0.,-1.],2:[0.,-2.,0.,2.],3:[-2.,0.,8.,0.,-6.],4:[0.,16.,0.,-40.,0.,24.]}
        ranges={}
        for j in range(1,order+3):
            pl,pu=polynomial_range(lo,hi,polys[j]);coef=a*w**j
            ranges[j]=(np.minimum(pl*coef,pu*coef).sum(axis=1),np.maximum(pl*coef,pu*coef).sum(axis=1))
        qr=q_range(left,right);qpr=(1-2*right,1-2*left);omega=self.frequency*np.pi
        if order==1:
            t3=interval_product(qr,ranges[3]);t2=interval_product(qpr,ranges[2])
            vl=t3[0]+3*t2[0]-6*ranges[1][1];vu=t3[1]+3*t2[1]-6*ranges[1][0]
            fl,fu=trig_range(left,right,omega,'cos');fl*=omega**3;fu*=omega**3
            rl,ru=-vu-fu,-vl-fl
        else:
            t4=interval_product(qr,ranges[4]);t3=interval_product(qpr,ranges[3])
            vl=t4[0]+4*t3[0]-12*ranges[2][1];vu=t4[1]+4*t3[1]-12*ranges[2][0]
            fl,fu=trig_range(left,right,omega,'sin');fl*=omega**4;fu*=omega**4
            rl,ru=-vu+fl,-vl+fu
        return np.maximum(abs(rl),abs(ru))
