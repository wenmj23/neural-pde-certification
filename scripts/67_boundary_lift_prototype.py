"""Research-only known nonzero boundary error prototype, no public CLI edits.
Boundary error is exact rational. Residual adapter represents the zero-boundary
part; affine lift is added to reconstructed error before norm enclosure.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,time
from math import comb
from flint import arb
from arb_reconstruction import norm_enclosure,point
from arb_shaped_remainder import ShapedVerifier,power_to_bernstein
from arb_cubic_reconstruction import CubicVerifier
from arb_certification import finite,exact_fraction
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261030/boundary_lift';OUT.mkdir(parents=True,exist_ok=False)
source=(ROOT/'scripts/61_exact_polynomial_audit.py').read_text().split('meshes=')[0];source=source.replace("OUT=Path(__file__).resolve().parents[2]/'run_20261030/polynomial_audit';OUT.mkdir(parents=True,exist_ok=False)","OUT=None");ns={'__file__':str(ROOT/'scripts/61_exact_polynomial_audit.py')};exec(compile(source,'polynomial_definitions','exec'),ns)
POLYS=ns['POLYS'];Polynomial=ns['Polynomial'];error_power=ns['error_power']
def reference(coeff,g0,g1,depth=11):
 p=error_power(coeff);p[0]+=g0;p[1]+=g1-g0;n=len(p)-1
 bern=[sum((p[j]*F(comb(k,j),comb(n,j)) for j in range(k+1)),F(0)) for k in range(n+1)];pieces=[bern]
 for _ in range(depth):
  out=[]
  for row in pieces:
   row=list(row);left=[row[0]];right=[row[-1]]
   while len(row)>1:
    row=[(a+b)/2 for a,b in zip(row[:-1],row[1:])];left.append(row[0]);right.append(row[-1])
   out.extend([left,list(reversed(right))])
  pieces=out
 return max(max(abs(p[0]),abs(p[-1])) for p in pieces),max(abs(c) for p in pieces for c in p)
def lift_norm(verifier,cells,left,right):
 # Rebuild the same signed interpolant z, adding the known affine lift.
 slope=point(right-left)
 for c in cells:
  a,b,h=c['a'],c['b'],c['h']
  if verifier.degree==2:
   r0,r1=c['r0'],c['r1'];mid=(a+b)/2
   term=h*((1-a)*r0+4*(1-mid)*c['rm']+(1-b)*r1)/6
  else:
   p=c['p'];mass=sum((v/(j+1) for j,v in enumerate(p)),arb(0));first=sum((v/(j+2) for j,v in enumerate(p)),arb(0));term=h*((1-a)*mass-h*first)
  slope=finite(slope+term,'lifted initial slope')
 value=point(left);lo=hi=arb(0)
 for c in cells:
  h=c['h']
  if verifier.degree==2:
   r0,r1=c['r0'],c['r1'];m=c['rm'];B=4*m-3*r0-r1;C=2*(r0+r1-2*m);h2=h*h
   coeff=[value,value+slope*h/4,value+slope*h/2-r0*h2/12,value+3*slope*h/4-r0*h2/4-B*h2/24,value+slope*h-r0*h2/2-B*h2/6-C*h2/12];mass=h*(r0+4*m+r1)/6
  else:
   p=c['p'];power=[value,slope*h]+[-h*h*v/((j+1)*(j+2)) for j,v in enumerate(p)];coeff=power_to_bernstein(power);mass=h*sum((v/(j+1) for j,v in enumerate(p)),arb(0))
  low,high=norm_enclosure(coeff,verifier.depth);lo=max(lo,low);hi=max(hi,high);value=coeff[-1];slope=finite(slope-mass,'lifted slope')
 if not value.contains(point(right)):raise ArithmeticError('Lifted endpoint misses boundary')
 return lo,hi
def certify(coeff,nodes,degree,g0,g1):
 problem=Polynomial(coeff);v=ShapedVerifier(problem,'moment') if degree==2 else CubicVerifier(problem,'moment');cells=v.prepare(nodes)
 # The original residual remainder bound applies unchanged because affine lift has zero second derivative.
 zl,zu=lift_norm(v,cells,g0,g1);_,_,du=v._norms(cells);lower=exact_fraction(max(arb(0),(zl-du).lower()));upper=exact_fraction((zu+du).upper())
 base=v.certify(nodes);L=exact_fraction(base['lower']);U=exact_fraction(base['upper']);A=max(abs(g0),abs(g1));triangle=(max(F(0),L-A,A-U),A+U)
 return lower,upper,exact_fraction(du),triangle
boundary={'zero':(F(0),F(0)),'same_positive':(F(1,10000),F(1,10000)),'opposite':(F(-1,10000),F(1,10000)),'unequal':(F(1,20000),F(-1,40000))}
meshes={'one':[F(0),F(1)],'four':[F(i,4) for i in range(5)],'sixteen':[F(i,16) for i in range(17)],'nonuniform':[F(0),F(1,17),F(1,3),F(7,10),F(1)]};rows=[]
for name,coeff in POLYS.items():
 for bname,(g0,g1) in boundary.items():
  rlo,rhi=reference(coeff,g0,g1)
  for mesh,nodes in meshes.items():
   for degree in (2,3):
    lo,hi,du,(tlo,thi)=certify(coeff,nodes,degree,g0,g1);compatible=lo<=rhi and hi>=rlo;overlap=max(lo,tlo)<=min(hi,thi)
    known=(max(abs(g0),abs(g1)) if name=='zero' else None);knownok=(lo<=known<=hi) if known is not None else None
    rows.append({'case':name,'boundary':bname,'g0':str(g0),'g1':str(g1),'mesh':mesh,'degree':degree,'lower':str(lo),'upper':str(hi),'remainder_upper':str(du),'reference_lower':str(rlo),'reference_upper':str(rhi),'compatible':compatible,'known_zero_residual_exact_contained':knownok,'triangle_lower':str(tlo),'triangle_upper':str(thi),'triangle_overlap':overlap,'lift_tighter_width':hi-lo<thi-tlo});(OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n')
 print(name,'done',flush=True)
summary={'synthetic_residuals':len(POLYS),'boundary_pairs':len(boundary),'cases':len(rows),'compatible':sum(x['compatible'] for x in rows),'triangle_overlap':sum(x['triangle_overlap'] for x in rows),'known_exact_checks':sum(x['known_zero_residual_exact_contained'] is not None for x in rows),'known_exact_passed':sum(x['known_zero_residual_exact_contained'] is True for x in rows),'lift_tighter_width':sum(x['lift_tighter_width'] for x in rows),'failures':[x for x in rows if not x['compatible'] or not x['triangle_overlap'] or x['known_zero_residual_exact_contained'] is False],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k!='failures'}),'failures',len(summary['failures']),flush=True)
