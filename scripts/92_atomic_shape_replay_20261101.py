"""Compare remainder propagation methods while reusing the same cell data."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import json,hashlib,os

def durable_json(path,value):
 data=(json.dumps(value,indent=2)+"\n").encode()
 temp=path.with_name(path.name+".tmp-"+str(os.getpid()))
 with temp.open("xb") as f:
  f.write(data);f.flush();os.fsync(f.fileno())
 os.replace(temp,path)

from arb_cubic_reconstruction import CubicVerifier
from arb_shaped_remainder import ShapedVerifier
from arb_certification import exact_fraction
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261101/shape_ablation_atomic';OUT.mkdir(parents=True,exist_ok=False)
s=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];s=s.replace("OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)","OUT=None");ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(s,'poly_forcing_definitions','exec'),ns)
fixture=ROOT.parent/'run_20261101/poly_new_training/frozen_polynomial_pinns.json';data=fixture.read_bytes();sha=hashlib.sha256(data).hexdigest();assert sha==(fixture.parent/'fixture_sha256.txt').read_text().strip();nets=json.loads(data)['networks'];refs=json.loads((ROOT.parent/'run_20261101/poly_new_certification/references.json').read_text());ablation=json.loads((ROOT.parent/'run_20261101/matched_residual_ablation/records.json').read_text());lookup={(x['network'],x['mesh_intervals'],x['method']):x for x in ablation};fr=ns['fr'];rows=[]
config={'fixture_sha256':sha,'grids':[16,32,64,128],'modes':{'quadratic':['constant','full','moment'],'cubic':['constant','moment']},'policy':'same problem, same mesh and cached source data within each degree; no timing comparison','reference_use':'validation only','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};durable_json(OUT/'config.json',config)
for net in nets:
 name=net['label'];coeff=net['source_coefficients'];ref=refs[name];rl,ru=fr(ref['lower']),fr(ref['upper'])
 for n in config['grids']:
  for deg,modes in ((2,('constant','full','moment')),(3,('constant','moment'))):
   p=ns['PolynomialForcing'](fixture,name,sha,coeff);v=ShapedVerifier(p,'moment') if deg==2 else CubicVerifier(p,'moment');v.uniform(n);v.prepare(v.nodes);counts=(v.residual_evaluations,v.curvature_evaluations)
   for mode in modes:
    v.mode=mode;r=v.certify();lo,hi,du=exact_fraction(r['lower']),exact_fraction(r['upper']),exact_fraction(r['remainder_upper']);old=lookup[(name,n,'signed_quadratic' if deg==2 else 'signed_cubic')]
    same=(lo==F(old['lower']) and hi==F(old['upper'])) if mode=='moment' else None
    row={'network':name,'source':net['source_name'],'depth':net['depth'],'seed':net['seed'],'intervals':n,'degree':deg,'mode':mode,'lower':str(lo),'upper':str(hi),'D_upper':str(du),'reference_overlap':lo<=ru and hi>=rl,'reference_contained':lo<=rl and hi>=ru,'same_as_previous_moment':same,'no_extra_model_calls':counts==(v.residual_evaluations,v.curvature_evaluations),'gap_1_02':(lo==hi==0) or (lo>0 and hi<=F(51,50)*lo)};rows.append(row)
  durable_json(OUT/'records.json',rows)
 print(name,'complete',flush=True)
comparisons=[]
for net in nets:
 for n in config['grids']:
  for deg in (2,3):
   z={x['mode']:x for x in rows if x['network']==net['label'] and x['intervals']==n and x['degree']==deg};a=F(z['constant']['D_upper']);b=F(z['moment']['D_upper']);c=F(z['full']['D_upper']) if deg==2 else None
   comparisons.append({'network':net['label'],'intervals':n,'degree':deg,'moment_less':b<a,'moment_equal':b==a,'moment_greater':b>a,'full_less_than_moment':c<b if c is not None else None,'full_greater_than_moment':c>b if c is not None else None,'moment_to_constant':float(b/a) if a else None})
summary={'networks':len(nets),'bounds':len(rows),'groups':len(comparisons),'reference_overlap':sum(x['reference_overlap'] for x in rows),'reference_contained':sum(x['reference_contained'] for x in rows),'cache_no_extra':sum(x['no_extra_model_calls'] for x in rows),'moment_exact_checks':sum(x['same_as_previous_moment'] is not None for x in rows),'moment_exact_matches':sum(x['same_as_previous_moment'] is True for x in rows),'by_degree':{str(d):{'comparisons':sum(x['degree']==d for x in comparisons),'moment_less':sum(x['degree']==d and x['moment_less'] for x in comparisons),'moment_equal':sum(x['degree']==d and x['moment_equal'] for x in comparisons),'moment_greater':sum(x['degree']==d and x['moment_greater'] for x in comparisons)} for d in (2,3)},'quadratic_full_tighter_than_moment':sum(x['full_less_than_moment'] is True for x in comparisons),'quadratic_full_looser_than_moment':sum(x['full_greater_than_moment'] is True for x in comparisons)}
durable_json(OUT/'comparisons.json',comparisons);durable_json(OUT/'summary.json',summary);print('SUMMARY',json.dumps(summary),flush=True)
if summary['reference_overlap']!=len(rows) or summary['cache_no_extra']!=len(rows) or summary['moment_exact_matches']!=summary['moment_exact_checks']:raise SystemExit(1)
