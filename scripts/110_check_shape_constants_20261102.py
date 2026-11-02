from fractions import Fraction as F
from math import comb

# Polynomial coefficients are in ascending powers, exact rational arithmetic.
def mul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c

def integ(p,a,b):
    return sum((v*(b**(i+1)-a**(i+1))/F(i+1) for i,v in enumerate(p)), F(0))

def poly(nodes):
    p=[F(1)]
    for n in nodes: p=mul(p,[-n,F(1)])
    return p

# Sign of nodal product alternates between successive roots.
def abs_integral(p,nodes):
    return sum((abs(integ(p,a,b)) for a,b in zip(nodes,nodes[1:])), F(0))

n2=[F(0),F(1,2),F(1)]
n3=[F(0),F(1,3),F(2,3),F(1)]
p2,p3=poly(n2),poly(n3)
a2,a3=abs_integral(p2,n2),abs_integral(p3,n3)
assert a2==F(1,32) and a3==F(49,7290)
assert 81*a3==F(49,90)
# Cubic maximum: in u=t-1/2, nodal product is
# (u²-1/4)(u²-1/36); its critical values are 1/144 and -1/81.
assert (F(0)-F(1,4))*(F(0)-F(1,36))==F(1,144)
assert (F(5,36)-F(1,4))*(F(5,36)-F(1,36))==-F(1,81)
# Independently integrate twice the quadratic interpolant r0+B*t+C*t²,
# then convert power coefficients to degree-four Bernstein coefficients.
r0,rm,r1=F(2,7),F(-3,8),F(5,11)
B=4*rm-3*r0-r1; C=2*(r0+r1-2*rm)
v,s,h=F(4,9),F(-1,5),F(3,7)
power=[v,s*h,-r0*h*h/2,-B*h*h/6,-C*h*h/12]
bern=[sum((power[j]*F(comb(k,j),comb(4,j)) for j in range(k+1)),F(0)) for k in range(5)]
expected=[v,v+s*h/4,v+s*h/2-r0*h*h/12,v+3*s*h/4-r0*h*h/4-B*h*h/24,v+s*h-r0*h*h/2-B*h*h/6-C*h*h/12]
assert bern==expected
print('quadratic |nodal product| integral:', a2, '(times 12 sqrt(3) gives 3 sqrt(3)/8)')
print('cubic |nodal product| integral:', a3, '(times 81 gives 49/90)')
print('cubic max nodal product: 1/81; interpolation factor 1/1944')
print('quadratic auxiliary Bernstein conversion: exact rational check passed')
