"""Exact-rational integration of interpolation-node shapes."""
from fractions import Fraction as F
from pathlib import Path
from math import comb
import hashlib,json
from arb_cubic_reconstruction import CUBIC_AVERAGE
OUT=Path(__file__).resolve().parents[2]/'run_20261030/exact_shape_identities';OUT.mkdir(parents=True,exist_ok=False)
def integ(p,a,b):return sum((c*(b**(i+1)-a**(i+1))/F(i+1) for i,c in enumerate(p)),F(0))
def shiftmul(p):return [F(0)]+p
def reflect(p):
 q=[F(0)]*len(p)
 for j,c in enumerate(p):
  for k in range(j+1):q[k]+=c*comb(j,k)*(-1)**k
 return q
q=[F(0),F(1,2),F(-3,2),F(1)]
c=[F(0),F(-2,9),F(11,9),F(-2),F(1)]
qpieces=[(F(0),F(1,2),1),(F(1,2),F(1),-1)]
cpieces=[(F(0),F(1,3),-1),(F(1,3),F(2,3),1),(F(2,3),F(1),-1)]
def abs_integral(p,pieces):return sum((s*integ(p,a,b) for a,b,s in pieces),F(0))
qmass=abs_integral(q,qpieces);qfirst=abs_integral(shiftmul(q),qpieces);cmass=81*abs_integral(c,cpieces);cfirst=81*abs_integral(shiftmul(c),cpieces)
# c(t)=(y-1/4)(y-1/36), y=(t-1/2)^2 in [0,1/4].
# Vertex y=5/36 yields -1/81; endpoints yield 1/144 and 0.
y=F(5,36);vertex=(y-F(1,4))*(y-F(1,36))
checks={'quadratic_symmetry':reflect(q)==[-x for x in q],'cubic_symmetry':reflect(c)==c,'quadratic_mass_1_32':qmass==F(1,32),'quadratic_first_half_mass':qfirst==qmass/2,'cubic_average_49_90':cmass==F(49,90),'cubic_first_half_mass':cfirst==cmass/2,'cubic_code_average_matches':cmass==CUBIC_AVERAGE,'cubic_vertex_minus_1_81':vertex==-F(1,81),'cubic_center_1_144':(F(-1,4)*F(-1,36))==F(1,144)}
summary={'checks':checks,'passed':sum(checks.values()),'total':len(checks),'quadratic_mass_without_12sqrt3':str(qmass),'cubic_mass':str(cmass),'cubic_first_moment':str(cfirst),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
