"""Prespecified uniform-grid certification; independent reference afterwards."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import json,hashlib,time,traceback
from arb_cubic_reconstruction import CubicVerifier
from arb_shaped_remainder import ShapedVerifier
from arb_certification import exact_fraction
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261101'/'poly_new_certification';OUT.mkdir(parents=True,exist_ok=False)
s=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];s=s.replace("OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)","OUT=None");ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(s,'poly_forcing_definitions','exec'),ns)
fixture=ROOT.parent/'run_20261101/poly_new_training/frozen_polynomial_pinns.json';data=fixture.read_bytes();sha=hashlib.sha256(data).hexdigest();assert sha==(fixture.parent/'fixture_sha256.txt').read_text().strip();nets=json.loads(data)['networks'];rows=[];refs={};fr=ns['fr']
config={'fixture_sha256':sha,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'prec_bits':103,'degrees':[2,3],'uniform_mesh_sequence':[16,32,64,128,256,512],'goal':'U <= 51/50 L with L>0 or exact zero','decision_tol':'1/1000','true_solution':'not used during certificate attempts, choices, or stopping; reference only after both degree runs per network','reference_target':'relative 1e-7, absolute 1e-14, cap 5000','baseline_commit':'8e09c023419ef92dfff9b12f29867c915ae0a4fd'};(OUT/'locked_config.json').write_text(json.dumps(config,indent=2)+'\n')
for net in nets:
 name=net['label'];pair=[]
 for degree in (2,3):
  tick=time.perf_counter();trace=[];final=None;status='budget_exhausted';err=None
  try:
   p=ns['PolynomialForcing'](fixture,name,sha,net['source_coefficients']);v=ShapedVerifier(p,'moment') if degree==2 else CubicVerifier(p,'moment')
   for n in config['uniform_mesh_sequence']:
    v.uniform(n);r=v.certify();lo,hi=exact_fraction(r['lower']),exact_fraction(r['upper']);met=bool((lo==0 and hi==0) or (lo>0 and hi<=F(51,50)*lo));entry={'intervals':n,'lower':str(lo),'upper':str(hi),'remainder_upper':str(exact_fraction(r['remainder_upper'])),'gap_1_02':met};trace.append(entry);final=entry
    if met:status='goal_met';break
  except Exception as e:status='exception';err={'message':str(e),'traceback':traceback.format_exc()}
  row={'network':name,'source':net['source_name'],'seed':net['seed'],'depth':net['depth'],'degree':degree,'status':status,'final':final,'trace':trace,'elapsed_seconds':time.perf_counter()-tick,'error':err};rows.append(row);pair.append(row);(OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n');print(name,'degree',degree,status,'intervals',final['intervals'] if final else None,flush=True)
 try:
  ref=ns['PolynomialReference'](net,net['source_coefficients']).maximum(relative=F(1,10**7),absolute=F(1,10**14),cap=5000);refs[name]=ref;rl,ru=fr(ref['lower']),fr(ref['upper'])
  for row in pair:
   f=row['final'];row['reference_converged']=ref['converged'];row['reference_overlap']=bool(f and F(f['lower'])<=ru and F(f['upper'])>=rl);row['reference_contained']=bool(f and F(f['lower'])<=rl and F(f['upper'])>=ru)
 except Exception as e:
  refs[name]={'converged':False,'error':str(e),'traceback':traceback.format_exc()}
  for row in pair:row['reference_converged']=False;row['reference_overlap']=None;row['reference_contained']=None
 (OUT/'references.json').write_text(json.dumps(refs,indent=2)+'\n');(OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n')
summary={'networks':len(nets),'requests':len(rows),'status_counts':dict(Counter(r['status'] for r in rows)),'reference_cases':len(refs),'reference_converged':sum(r.get('converged',False) for r in refs.values()),'reference_overlap':sum(r.get('reference_overlap') is True for r in rows),'reference_contained':sum(r.get('reference_contained') is True for r in rows),'noncontained':[{'network':r['network'],'degree':r['degree'],'overlap':r.get('reference_overlap')} for r in rows if r.get('reference_contained') is False],'methods':{str(d):{'requests':sum(r['degree']==d for r in rows),'goal_met':sum(r['degree']==d and r['status']=='goal_met' for r in rows)} for d in (2,3)}};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps(summary),flush=True)
