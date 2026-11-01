"""Validate upper bounds and account for cumulative derivative-work counters."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent/'run_20261101';OUT=BASE/'benchmark_trace_checks_corrected';OUT.mkdir(parents=True,exist_ok=False)
a=json.loads((BASE/'matched_residual_ablation/records.json').read_text());lookup={(x['network'],x['method'],x['mesh_intervals']):x for x in a};b=json.loads((BASE/'matched_acceptance_benchmark/records.json').read_text());checks=[]
for r in b:
 cumulative=0
 for t in r['trace']:
  k=(r['network'],r['method'],t['intervals']);old=lookup[k];cumulative+=old['derivative_evaluations'];checks.append({'network':r['network'],'method':r['method'],'tolerance':r['tolerance'],'repeat':r['repeat'],'intervals':t['intervals'],'upper_exact_match':F(t['upper'])==F(old['upper']),'residual_calls_match':t['residual_evaluations']==old['residual_evaluations'],'derivative_calls_match_cumulative':t['derivative_evaluations']==cumulative,'fresh_single_mesh_derivative_calls':old['derivative_evaluations'],'cumulative_expected':cumulative,'cumulative_observed':t['derivative_evaluations']})
summary={'request_count':len(b),'trace_checks':len(checks),'upper_matches':sum(x['upper_exact_match'] for x in checks),'residual_call_matches':sum(x['residual_calls_match'] for x in checks),'derivative_cumulative_matches':sum(x['derivative_calls_match_cumulative'] for x in checks),'failures':[x for x in checks if not(x['upper_exact_match'] and x['residual_calls_match'] and x['derivative_calls_match_cumulative'])],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'prior_check_issue':'script94 wrongly compared cumulative derivative count to fresh single-mesh count; saved prior results unchanged'};(OUT/'records.json').write_text(json.dumps(checks,indent=2)+'\n');(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='failures'}),'failures',len(summary['failures']));raise SystemExit(1 if summary['failures'] else 0)
