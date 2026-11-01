"""Matched quadratic-mode timing: same frozen network, mesh and source."""
from pathlib import Path
from fractions import Fraction as F
from statistics import median
import json,hashlib,time
from arb_shaped_remainder import ShapedVerifier
from arb_certification import exact_fraction
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261101/shaped_remainder_timing';OUT.mkdir(parents=True,exist_ok=False)
s=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];s=s.replace("OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)","OUT=None");ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(s,'poly_forcing_definitions','exec'),ns)
fixture=ROOT.parent/'run_20261101/poly_new_training/frozen_polynomial_pinns.json';data=fixture.read_bytes();sha=hashlib.sha256(data).hexdigest();nets=json.loads(data)['networks'];selected=[n for n in nets if n['seed']==203];prior=json.loads((ROOT.parent/'run_20261101/shaped_remainder_ablation/records.json').read_text());lookup={(x['network'],x['intervals'],x['degree'],x['mode']):x for x in prior};rows=[]
for net in selected:
 for n in (32,64):
  for rep in range(1,4):
   modes=('constant','full','moment') if rep%2 else ('moment','full','constant')
   for mode in modes:
    start=time.perf_counter();p=ns['PolynomialForcing'](fixture,net['label'],sha,net['source_coefficients']);v=ShapedVerifier(p,mode);v.uniform(n);r=v.certify();lo,hi,du=exact_fraction(r['lower']),exact_fraction(r['upper']),exact_fraction(r['remainder_upper']);elapsed=time.perf_counter()-start;old=lookup[(net['label'],n,2,mode)];ok=lo==F(old['lower']) and hi==F(old['upper']) and du==F(old['D_upper']);rows.append({'network':net['label'],'intervals':n,'mode':mode,'repeat':rep,'seconds':elapsed,'exact_match':ok});(OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n')
  print(net['label'],n,'done',flush=True)
groups=[]
for net in selected:
 for n in (32,64):
  med={m:median(x['seconds'] for x in rows if x['network']==net['label'] and x['intervals']==n and x['mode']==m) for m in ('constant','full','moment')};groups.append({'network':net['label'],'intervals':n,'median_seconds':med,'full_vs_constant':med['full']/med['constant'],'moment_vs_constant':med['moment']/med['constant'],'full_vs_moment':med['full']/med['moment']})
summary={'networks':len(selected),'groups':len(groups),'runs':len(rows),'exact_matches':sum(x['exact_match'] for x in rows),'median_of_group_full_vs_constant':median(x['full_vs_constant'] for x in groups),'median_of_group_moment_vs_constant':median(x['moment_vs_constant'] for x in groups),'median_of_group_full_vs_moment':median(x['full_vs_moment'] for x in groups),'groups_detail':groups,'scope':'single process and machine; per run includes fixture/model construction, derivatives, certification; excludes imports, output and reference; three repeats','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k!='groups_detail'}),flush=True)
if summary['exact_matches']!=len(rows):raise SystemExit(1)
