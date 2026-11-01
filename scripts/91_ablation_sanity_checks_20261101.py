"""Known polynomial problems sanity-check shared-data residual ablation."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
from arb_certification import exact_fraction
from arb_shaped_remainder import ShapedVerifier
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261101/ablation_sanity';OUT.mkdir(parents=True,exist_ok=False)
s=(ROOT/'scripts/61_exact_polynomial_audit.py').read_text().split('meshes=')[0];s=s.replace("OUT=Path(__file__).resolve().parents[2]/'run_20261030/polynomial_audit';OUT.mkdir(parents=True,exist_ok=False)","OUT=None");n={'__file__':str(ROOT/'scripts/61_exact_polynomial_audit.py')};exec(compile(s,'poly_definitions','exec'),n)
s=(ROOT/'scripts/88_matched_residual_ablation_20261101.py').read_text().split('config=')[0];s=s.replace("OUT=ROOT.parent/'run_20261101/matched_residual_ablation';OUT.mkdir(parents=True,exist_ok=False)","OUT=None");a={'__file__':str(ROOT/'scripts/88_matched_residual_ablation_20261101.py')};exec(compile(s,'ablation_definitions','exec'),a);rows=[]
for name in ('zero','constant','negative_constant','node_hidden_quartic'):
 rl,ru=n['reference'](n['POLYS'][name],depth=12)
 for intervals in (1,4,16,32):
  p=n['Polynomial'](n['POLYS'][name]);v=ShapedVerifier(p,'moment');v.uniform(intervals);r=v.certify();calls=(v.residual_evaluations,v.curvature_evaluations);glob,local=a['bounds'](v.prepare(v.nodes),v.depth);sg=exact_fraction(r['upper']);correct=all(u>=rl for u in (glob,local,sg));exact=(F(0) if name=='zero' else F(1,80000) if name in ('constant','negative_constant') else None);saturation=all(0<=u-exact<F(1,10**25) for u in (glob,local,sg)) if exact is not None else None
  rows.append({'case':name,'intervals':intervals,'global_upper':str(glob),'spatial_upper':str(local),'signed_upper':str(sg),'reference_lower':str(rl),'reference_upper':str(ru),'covers_reference_lower':correct,'exact_saturation':saturation,'no_extra_model_calls':calls==(v.residual_evaluations,v.curvature_evaluations)})
summary={'groups':len(rows),'coverage_checks':sum(x['covers_reference_lower'] for x in rows),'exact_saturation_checks':sum(x['exact_saturation'] is not None for x in rows),'exact_saturation_passed':sum(x['exact_saturation'] is True for x in rows),'cache_checks_passed':sum(x['no_extra_model_calls'] for x in rows),'failures':[x for x in rows if not x['covers_reference_lower'] or x['exact_saturation'] is False or not x['no_extra_model_calls']],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(OUT/'records.json').write_text(json.dumps(rows,indent=2)+'\n');(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2));raise SystemExit(1 if summary['failures'] else 0)
