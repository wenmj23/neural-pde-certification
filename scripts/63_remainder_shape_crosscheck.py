"""Numerical Green integration of reconstructed remainder shapes; diagnostic only."""
import importlib.util,json,hashlib
from pathlib import Path
from fractions import Fraction as F
import mpmath as mp
from arb_shaped_remainder import ShapedVerifier
from arb_cubic_reconstruction import CubicVerifier
from arb_certification import exact_fraction
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261030/remainder_shape_check';OUT.mkdir(parents=True,exist_ok=False)
# Loading definitions without rerunning script 61's driver.
source=(ROOT/'scripts/61_exact_polynomial_audit.py').read_text();prefix=source.split('meshes=')[0];prefix=prefix.replace("OUT=Path(__file__).resolve().parents[2]/'run_20261030/polynomial_audit';OUT.mkdir(parents=True,exist_ok=False)","OUT=None")
ns={'__file__':str(ROOT/'scripts/61_exact_polynomial_audit.py')};exec(compile(prefix,'polynomial_definitions','exec'),ns)
mp.mp.dps=70
meshes={'one':[F(0),F(1)],'four':[F(i,4) for i in range(5)],'nonuniform':[F(0),F(1,17),F(1,3),F(7,10),F(1)]}
names=['constant','linear','quadratic','cubic','node_hidden_quartic','alternating_sextic','asymmetric_sextic']; rows=[]
def mpq(x):return mp.mpf(x.numerator)/x.denominator
def compute_max(cells,nodes,degree):
 pieces=[]
 for c,a,b in zip(cells,nodes[:-1],nodes[1:]):
  aa=mpq(a);h=mpq(b-a);delta=mpq(exact_fraction(c['delta']))
  if degree==2: roots=[mp.mpf(0),mp.mpf(1)/2,mp.mpf(1)]; signs=[1,-1];co=[0,mp.mpf(1)/2,-mp.mpf(3)/2,1];amp=12*mp.sqrt(3)*delta
  else: roots=[mp.mpf(0),mp.mpf(1)/3,mp.mpf(2)/3,mp.mpf(1)];signs=[-1,1,-1];co=[0,-mp.mpf(2)/9,mp.mpf(11)/9,-2,1];amp=81*delta
  for l,r,s in zip(roots[:-1],roots[1:],signs):pieces.append((aa,h,l,r,s,co,amp))
 def integ(piece,limit):
  a,h,l,r,s,co,amp=piece;t=min(r,max(l,(limit-a)/h))
  if t<=l:return mp.mpf(0),mp.mpf(0)
  mass=amp*s*h*sum(co[j]*(t**(j+1)-l**(j+1))/(j+1) for j in range(len(co)))
  first=amp*s*h*sum(co[j]*(a*(t**(j+1)-l**(j+1))/(j+1)+h*(t**(j+2)-l**(j+2))/(j+2)) for j in range(len(co)))
  return mass,first
 mass0=first0=mp.mpf(0)
 for p in pieces:
  m,f=integ(p,mp.mpf(1));mass0+=m;first0+=f
 slope=mass0-first0
 if mass0==0:return mp.mpf(0)
 l=mp.mpf(0);r=mp.mpf(1)
 for _ in range(140):
  x=(l+r)/2;m=sum((integ(p,x)[0] for p in pieces),mp.mpf(0))
  if m<slope:l=x
  else:r=x
 x=(l+r)/2;m=f=mp.mpf(0)
 for p in pieces:
  a,b=integ(p,x);m+=a;f+=b
 return x*slope-x*m+f
for name in names:
 for mesh,nodes in meshes.items():
  for degree,mode in [(2,'constant'),(2,'full'),(2,'moment'),(3,'constant'),(3,'moment')]:
   prob=ns['Polynomial'](ns['POLYS'][name]);v=ShapedVerifier(prob,mode) if degree==2 else CubicVerifier(prob,mode);cells=v.prepare(nodes);numerical=compute_max(cells,nodes,degree);r=v.certify(nodes);upper=exact_fraction(r['remainder_upper']);u=mpq(upper)
   # Numerical comparisons allow only a tiny diagnostic tolerance, not a rigorous rounding guarantee.
   ok=numerical <= u+mp.mpf('1e-58')
   rows.append({'case':name,'mesh':mesh,'degree':degree,'mode':mode,'green_max_numeric':mp.nstr(numerical,60),'reported_D_upper':str(upper),'diagnostic_passed':bool(ok)});(OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n')
 print(name,'done',flush=True)
summary={'configurations':len(rows),'passed':sum(x['diagnostic_passed'] for x in rows),'failures':[x for x in rows if not x['diagnostic_passed']],'caveat':'mpmath approximate integration and bisection, not formal proof','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k!='failures'}),'failures',len(summary['failures']),flush=True)
