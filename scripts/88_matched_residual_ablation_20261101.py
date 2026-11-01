"""Matched quadratic-data ablation on Nov 1 networks; no runtime claim."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import json,hashlib
from flint import arb
from arb_certification import finite,exact_fraction
from arb_reconstruction import norm_enclosure
from arb_shaped_remainder import ShapedVerifier
from arb_cubic_reconstruction import CubicVerifier
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261101/matched_residual_ablation';OUT.mkdir(parents=True,exist_ok=False)
s=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];s=s.replace("OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)","OUT=None");ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(s,'poly_forcing_definitions','exec'),ns)
fixture=ROOT.parent/'run_20261101/poly_new_training/frozen_polynomial_pinns.json';data=fixture.read_bytes();sha=hashlib.sha256(data).hexdigest();assert sha==(fixture.parent/'fixture_sha256.txt').read_text().strip();nets=json.loads(data)['networks'];refs=json.loads((ROOT.parent/'run_20261101/poly_new_certification/references.json').read_text());fr=ns['fr'];rows=[]
def bounds(cells,depth):
 envelopes=[]
 for c in cells:
  r0,rm,r1=c['r0'],c['rm'],c['r1'];coef=[r0,2*rm-(r0+r1)/2,r1];_,pu=norm_enclosure(coef,depth);envelopes.append(finite(pu+c['delta'],'residual envelope').upper())
 glob=finite(max(envelopes)/8,'global residual bound').upper();slope=arb(0)
 for c,rho in zip(cells,envelopes):slope=finite(slope+rho*c['h']*(1-(c['a']+c['b'])/2),'initial slope')
 val=upper=arb(0)
 for c,rho in zip(cells,envelopes):
  h=c['h'];poly=[val,val+slope*h/2,val+slope*h-rho*h*h/2];_,u=norm_enclosure(poly,depth);upper=max(upper,u);val=poly[-1];slope=finite(slope-rho*h,'slope update')
 if not val.contains(0):raise ArithmeticError('Envelope endpoint misses zero')
 return exact_fraction(glob),exact_fraction(finite(upper,'spatial bound').upper())
def classification(lo,hi,tau):return 'within' if hi<=tau else 'exceeds' if lo>tau else 'unknown'
config={'fixture_sha256':sha,'networks':len(nets),'grid_sizes':[16,32,64,128],'quadratic_subdivision':2,'cubic_subdivision':2,'tolerances':['1/100000','1/1000000'],'scope':'4 upper bounds at same mesh; first 3 share quadratic evaluations; cubic uses more nodes and different derivatives; no speed comparison','reference_use':'validation only; no selection or stopping','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(OUT/'config.json').write_text(json.dumps(config,indent=2)+'\n')
for net in nets:
 name=net['label'];coeff=net['source_coefficients'];ref=refs[name];rl,ru=fr(ref['lower']),fr(ref['upper'])
 for n in config['grid_sizes']:
  p=ns['PolynomialForcing'](fixture,name,sha,coeff);v=ShapedVerifier(p,'moment');v.uniform(n);r=v.certify();counts=(v.residual_evaluations,v.curvature_evaluations);cells=v.prepare(v.nodes);glob,spatial=bounds(cells,v.depth);cache_unchanged=counts==(v.residual_evaluations,v.curvature_evaluations);ql,qu=exact_fraction(r['lower']),exact_fraction(r['upper'])
  p3=ns['PolynomialForcing'](fixture,name,sha,coeff);v3=CubicVerifier(p3,'moment');v3.uniform(n);r3=v3.certify();cl,cu=exact_fraction(r3['lower']),exact_fraction(r3['upper'])
  M={'global':(F(0),glob),'spatial':(F(0),spatial),'signed_quadratic':(ql,qu),'signed_cubic':(cl,cu)}
  for method,(lo,hi) in M.items():
   midpoint=(rl+ru)/2
   record={'network':name,'source':net['source_name'],'depth':net['depth'],'seed':net['seed'],'mesh_intervals':n,'method':method,'lower':str(lo),'upper':str(hi),'reference_lower':str(rl),'reference_upper':str(ru),'reference_converged':ref['converged'],'reference_overlap':lo<=ru and hi>=rl,'reference_contained':lo<=rl and hi>=ru,'upper_over_reference_midpoint':float(hi/midpoint) if midpoint else None,'cache_unchanged':cache_unchanged if method in ('global','spatial','signed_quadratic') else None,'residual_evaluations':counts[0] if method!='signed_cubic' else v3.residual_evaluations,'derivative_evaluations':counts[1] if method!='signed_cubic' else v3.curvature_evaluations,'decisions':{str(t):classification(lo,hi,t) for t in (F(1,100000),F(1,1000000))}};rows.append(record)
  (OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n')
 print(name,'completed',flush=True)
summary={'networks':len(nets),'mesh_sizes':len(config['grid_sizes']),'groups':len(rows)//4,'bounds':len(rows),'reference_overlap':sum(x['reference_overlap'] for x in rows),'reference_contained':sum(x['reference_contained'] for x in rows),'shared_quadratic_data_checks':sum(x['cache_unchanged'] is True for x in rows),'incompatible':[{'network':x['network'],'method':x['method'],'mesh_intervals':x['mesh_intervals']} for x in rows if not x['reference_overlap'] or not x['reference_converged']],'by_method':{method:{'count':sum(x['method']==method for x in rows),'decisions':{tol:dict(Counter(x['decisions'][tol] for x in rows if x['method']==method)) for tol in ('1/100000','1/1000000')}} for method in ('global','spatial','signed_quadratic','signed_cubic')}};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps(summary),flush=True)
if summary['incompatible']:raise SystemExit(1)
