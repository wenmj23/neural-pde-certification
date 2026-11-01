"""Fixed thresholds, budget-limited decision requests on frozen new networks."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import json,hashlib,time,traceback
from arb_cubic_reconstruction import CubicVerifier
from arb_shaped_remainder import ShapedVerifier
from arb_certification import exact_fraction
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261101'/'poly_fixed_tolerance';OUT.mkdir(parents=True,exist_ok=False)
s=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];s=s.replace("OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)","OUT=None");ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(s,'poly_forcing_definitions','exec'),ns)
fixture=ROOT.parent/'run_20261101/poly_new_training/frozen_polynomial_pinns.json';data=fixture.read_bytes();sha=hashlib.sha256(data).hexdigest();assert sha==(fixture.parent/'fixture_sha256.txt').read_text().strip();nets=json.loads(data)['networks'];refs=json.loads((ROOT.parent/'run_20261101/poly_new_certification/references.json').read_text());fr=ns['fr'];rows=[]
tols=[F(1,1000),F(1,100000),F(1,1000000)];config={'fixture_sha256':sha,'degrees':[2,3],'tolerances':[str(x) for x in tols],'max_intervals':128,'uniform_mesh_sequence':[16,32,64,128],'reference_policy':'existing references are read after each decision; not used in stop conditions','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'note':'thresholds fixed for this analysis but references existed beforehand; not a formally preregistered blind test'};(OUT/'config.json').write_text(json.dumps(config,indent=2)+'\n')
for net in nets:
 for degree in (2,3):
  for tau in tols:
   tick=time.perf_counter();trace=[];decision='unknown';err=None
   try:
    p=ns['PolynomialForcing'](fixture,net['label'],sha,net['source_coefficients']);v=ShapedVerifier(p,'moment') if degree==2 else CubicVerifier(p,'moment')
    for n in config['uniform_mesh_sequence']:
     v.uniform(n);r=v.certify();lo,hi=exact_fraction(r['lower']),exact_fraction(r['upper']);decision='within' if hi<=tau else 'exceeds' if lo>tau else 'unknown';trace.append({'intervals':n,'lower':str(lo),'upper':str(hi),'decision':decision})
     if decision!='unknown':break
   except Exception as e:decision='exception';err={'message':str(e),'traceback':traceback.format_exc()}
   # This cross-check cannot select or modify the request and occurs only now.
   ref=refs[net['label']];rl,ru=fr(ref['lower']),fr(ref['upper']);refdec='within' if ru<=tau else 'exceeds' if rl>tau else 'reference_unknown';compatible=decision=='unknown' or (decision==refdec and ref.get('converged',False));row={'network':net['label'],'source':net['source_name'],'seed':net['seed'],'depth':net['depth'],'degree':degree,'tolerance':str(tau),'decision':decision,'reference_decision':refdec,'reference_converged':ref.get('converged',False),'compatible':compatible,'trace':trace,'elapsed_seconds':time.perf_counter()-tick,'error':err};rows.append(row);(OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n')
 print(net['label'],'completed',flush=True)
summary={'networks':len(nets),'requests':len(rows),'outcomes':dict(Counter(x['decision'] for x in rows)),'compatible':sum(x['compatible'] for x in rows),'failures':[{'network':x['network'],'degree':x['degree'],'tolerance':x['tolerance'],'decision':x['decision'],'reference':x['reference_decision']} for x in rows if not x['compatible']],'by_tolerance':{str(t):dict(Counter(x['decision'] for x in rows if x['tolerance']==str(t))) for t in tols},'by_degree':{str(d):dict(Counter(x['decision'] for x in rows if x['degree']==d)) for d in (2,3)}};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps(summary),flush=True)
