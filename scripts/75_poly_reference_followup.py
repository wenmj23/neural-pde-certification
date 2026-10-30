"""Post-hoc tighter independent reference for the non-contained brackets."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261030/poly_reference_followup';OUT.mkdir(parents=True,exist_ok=False)
source=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];source=source.replace("OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)","OUT=None");ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(source,'poly_forcing_definitions','exec'),ns)
name='depth1_k2_seed49';source_name='linear_6x';fixture=ROOT/'configs/frozen_networks_holdout_2026-10-30.json';net=next(x for x in json.loads(fixture.read_text())['networks'] if x['label']==name);coeff=ns['SOURCES'][source_name];ref=ns['PolynomialReference'](net,coeff).maximum(relative=F(1,10**12),absolute=F(1,10**16),cap=5000);prior=json.loads((ROOT.parent/'run_20261030/polynomial_forcing/records.json').read_text());oldrefs=json.loads((ROOT.parent/'run_20261030/polynomial_forcing/references.json').read_text());old=oldrefs[name+'__'+source_name];fr=ns['fr'];rl,ru=fr(ref['lower']),fr(ref['upper']);checks=[]
for r in prior:
 if r['network']==name and r['source']==source_name:
  lo,hi=F(r['lower']),F(r['upper']);checks.append({'degree':r['degree'],'intervals':r['intervals'],'overlap':lo<=ru and hi>=rl,'contains_tight_reference':lo<=rl and hi>=ru,'original_contains_reference':r['reference_contained']})
summary={'network':name,'source':source_name,'relative_target':'1e-12 (post-hoc)','original_reference_width':str(fr(old['upper'])-fr(old['lower'])),'new_reference_width':str(ru-rl),'new_reference':ref,'checks':checks,'passed':sum(c['contains_tight_reference'] for c in checks),'count':len(checks),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps(summary),flush=True)
