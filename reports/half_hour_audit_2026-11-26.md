# 2026-11-26 半点任务：尖峰导数公式审阅

研究对象为已知最大误差的人工高斯尖峰，不是训练网络。对高斯导数多项式 0 至 6 阶和乘积误差 2 至 6 阶，在符号变量下核验，12 个恒等式通过；故意修改六阶常数项的控制被识别。结果为符号代数核查，依赖 SymPy 和源码提取正确，并不证明区间后端或连续证书整体正确。

输入源文件来自当天整点任务源码快照并经其清单核验；本材料仅保存该项审阅所需文件，不代表同步了所有半点和整点实验。整点最近一轮有两条证书与旧参考部分重合而未完整包含，原因尚待更紧的参考核查。用户指定的整点成果应由整点任务日终补充；不将该计划写成完成。

复现：将以下两个源码代码块、audit.py 保存到同一新目录，以安装 SymPy 的 Python 运行 audit.py。它生成 results.json；原始结果代码块供比对。不要在包含历史文件的目录中运行，因为脚本会写 results.json。原始过程与独立读回在本地临时目录，跨运行持久性未保证。此前在线追踪清单属于另一源目录，已按配置指向核验，48 项一致，不能把它当本地目录清单。历史一次只读检查误重写旧 readback 及历史数据不完整根因未明的事实仍需保留。

## arb_problem_families.py

SHA-256 `ec72ac0d015bee596abd52b295d87f40af899cdb16e87573e0de872c31ed7fab`

```python
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

```

## arb_cubic_reconstruction.py

SHA-256 `8d842983e88cfa1b1c9e15c63d85bb527aa1d52d40de1ea576a312075d91026d`

```python
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

```

## audit.py

SHA-256 `3fc773abcd747c0d467ce04cfcdd8224f8408eaaafcd5f38a5b21ced1d7d4e3c`

```python
import sympy as s,json,ast
from pathlib import Path
p=Path(__file__).parent
# Extract literal coefficient arrays without importing original code.
t=ast.parse((p/'arb_problem_families.py').read_text());cubic=ast.parse((p/'arb_cubic_reconstruction.py').read_text());coef=None;extra=None
for n in ast.walk(t):
 if isinstance(n,ast.Assign) and any(isinstance(v,ast.Attribute) and v.attr=='coefficients' for v in n.targets):
  try:
   candidate=ast.literal_eval(n.value)
   if isinstance(candidate,list) and candidate and isinstance(candidate[0],list):coef=candidate
  except (ValueError,TypeError):pass
for n in ast.walk(cubic):
 if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='append' and n.args:
  try:
   q=ast.literal_eval(n.args[0])
   if q==[-120,0,720,0,-480,0,64]:extra=q
  except (ValueError,TypeError):pass
assert coef is not None and extra is not None
coef.append(extra);tvar,x,h=s.symbols('t x h',real=True,nonzero=True);r=s.Integer(1);rows=[];rs=[]
for j in range(7):
 src=sum(a*tvar**k for k,a in enumerate(coef[j]));rows.append({'kind':'gaussian_polynomial','order':j,'difference':str(s.expand(src-r)),'pass':s.expand(src-r)==0});rs.append(r);r=s.diff(r,tvar)-2*tvar*r
G=s.exp(-((x-s.Rational(1,2))/h)**2);q=x*(1-x)
for n in range(2,7):
 exact=s.diff(q*G,x,n)/G;u=(x-s.Rational(1,2))/h
 impl=q*rs[n].subs(tvar,u)/h**n+n*(1-2*x)*rs[n-1].subs(tvar,u)/h**(n-1)-n*(n-1)*rs[n-2].subs(tvar,u)/h**(n-2)
 diff=s.simplify(exact-impl);rows.append({'kind':'leibniz_assembly','order':n,'difference':str(diff),'pass':diff==0})
wrong=coef[6].copy();wrong[0]+=1;fault=s.expand(sum(a*tvar**k for k,a in enumerate(wrong))-rs[6]);rows.append({'kind':'negative_control','difference':str(fault),'pass':fault!=0})
(p/'results.json').write_text(json.dumps(rows,indent=2)+'\n');print('records',len(rows),'passed',sum(a['pass'] for a in rows),'fault_diff',fault)

```

## results.json

SHA-256 `4cb3909f5dae6fcb313afeda9f08f240aee482fb130dca5d54e4751c9db32fab`

```json
[
  {
    "kind": "gaussian_polynomial",
    "order": 0,
    "difference": "0",
    "pass": true
  },
  {
    "kind": "gaussian_polynomial",
    "order": 1,
    "difference": "0",
    "pass": true
  },
  {
    "kind": "gaussian_polynomial",
    "order": 2,
    "difference": "0",
    "pass": true
  },
  {
    "kind": "gaussian_polynomial",
    "order": 3,
    "difference": "0",
    "pass": true
  },
  {
    "kind": "gaussian_polynomial",
    "order": 4,
    "difference": "0",
    "pass": true
  },
  {
    "kind": "gaussian_polynomial",
    "order": 5,
    "difference": "0",
    "pass": true
  },
  {
    "kind": "gaussian_polynomial",
    "order": 6,
    "difference": "0",
    "pass": true
  },
  {
    "kind": "leibniz_assembly",
    "order": 2,
    "difference": "0",
    "pass": true
  },
  {
    "kind": "leibniz_assembly",
    "order": 3,
    "difference": "0",
    "pass": true
  },
  {
    "kind": "leibniz_assembly",
    "order": 4,
    "difference": "0",
    "pass": true
  },
  {
    "kind": "leibniz_assembly",
    "order": 5,
    "difference": "0",
    "pass": true
  },
  {
    "kind": "leibniz_assembly",
    "order": 6,
    "difference": "0",
    "pass": true
  },
  {
    "kind": "negative_control",
    "difference": "1",
    "pass": true
  }
]

```

