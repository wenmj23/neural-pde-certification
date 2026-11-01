"""Recompute selected fixed grids; compare all corresponding adaptive traces."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
from arb_cubic_reconstruction import CubicVerifier
from arb_shaped_remainder import ShapedVerifier
from arb_certification import exact_fraction
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent/'run_20261101';OUT=BASE/'holdout_fresh_grid_checks';OUT.mkdir(parents=True,exist_ok=False)
s=(ROOT/'scripts/97_holdout_matched_benchmark_20261101.py').read_text().split('config={')[0];s=s.replace("OUT=ROOT.parent/'run_20261101/holdout_benchmark';OUT.mkdir(parents=True,exist_ok=False)","OUT=None");env={'__file__':str(ROOT/'scripts/97_holdout_matched_benchmark_20261101.py')};exec(compile(s,'benchmark_definitions','exec'),env)
nets=[n for n in env['nets'] if n['seed']==305];fixture=env['fixture'];sha=env['sha'];methods=env['methods'];ns=env['ns'];rows=json.loads((BASE/'holdout_benchmark/records.json').read_text());fresh={};evaluations=[]
for net in nets:
 for method in methods:
  for n in (16,32,64,128):
   p=ns['PolynomialForcing'](fixture,net['label'],sha,net['source_coefficients']);v=CubicVerifier(p,'moment') if method=='signed_cubic' else ShapedVerifier(p,'moment');v.uniform(n)
   u=env['upper_for_cells'](v.prepare(v.nodes),v.depth,method) if method in ('global','spatial') else exact_fraction(v.certify()['upper']);key=(net['label'],method,n);fresh[key]={'upper':str(u),'residual_calls':v.residual_evaluations,'derivative_calls':v.curvature_evaluations};evaluations.append({'network':net['label'],'method':method,'intervals':n,**fresh[key]})
 print(net['label'],'recomputed',flush=True)
checks=[]
for r in rows:
 if not any(n['label']==r['network'] for n in nets):continue
 cumulative=0
 for t in r['trace']:
  item=fresh[(r['network'],r['method'],t['intervals'])];cumulative+=item['derivative_calls'];checks.append({'network':r['network'],'method':r['method'],'tolerance':r['tolerance'],'repeat':r['repeat'],'intervals':t['intervals'],'upper_match':F(t['upper'])==F(item['upper']),'residual_calls_match':t['residual_evaluations']==item['residual_calls'],'cumulative_derivatives_match':t['derivative_evaluations']==cumulative})
summary={'networks':len(nets),'fresh_fixed_mesh_evaluations':len(evaluations),'trace_checks':len(checks),'upper_matches':sum(c['upper_match'] for c in checks),'residual_matches':sum(c['residual_calls_match'] for c in checks),'derivative_cumulative_matches':sum(c['cumulative_derivatives_match'] for c in checks),'failures':[c for c in checks if not(c['upper_match'] and c['residual_calls_match'] and c['cumulative_derivatives_match'])],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(OUT/'evaluations.json').write_text(json.dumps(evaluations,indent=2)+'\n');(OUT/'checks.json').write_text(json.dumps(checks,indent=2)+'\n');(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k!='failures'}),'failures',len(summary['failures']),flush=True)
if summary['failures']:raise SystemExit(1)
