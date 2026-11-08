"""Read back the post-hoc sine extension and matched baseline report values."""
from fractions import Fraction as F
from statistics import median
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results/2026-11-07';D=B/'analytic_sine32_extended';BASE=ROOT/'results/2026-11-02/new_pinn_baselines/accuracy';OUT=B/'report_evidence_readback.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 if OUT.exists():raise FileExistsError(OUT)
 s=json.loads((D/'summary.json').read_text());rows=json.loads((D/'records.json').read_text())
 if sha(D/'records.json')!=s['records_sha256'] or sha(ROOT/'scripts/265_sine32_budget_followup_20261107.py')!=s['script_sha256'] or sha(B/'analytic_report_controls/records.json')!=s['prior_records_sha256']:raise ValueError('digest')
 seen=set();extension=[]
 for x in rows:
  key=x['mesh'],x['degree'];L,U,K=map(F,(x['lower'],x['upper'],x['known_error']))
  if key in seen or K!=F(1,10000) or not 0<L<=K<=U or not x['contains'] or F(x['ratio'])!=U/L or x['ratio_met']!=(U<=F(51,50)*L):raise ValueError('extension arithmetic')
  seen.add(key);extension.append({'mesh':key[0],'degree':key[1],'ratio':float(U/L),'goal_met':x['ratio_met']})
 if seen!={(m,d) for m in (512,1024) for d in (2,3)} or s['records']!=4 or s['contained']!=4:raise ValueError('extension manifest')
 bs=json.loads((BASE/'summary.json').read_text());br=json.loads((BASE/'records.json').read_text())
 if sha(BASE/'records.json')!=bs['records_sha256'] or len(br)!=144:raise ValueError('baseline digest')
 comparisons=[]
 for method in ('global','spatial','signed_quadratic','signed_cubic'):
  group=[r for r in br if r['intervals']==64 and r['method']==method]
  if len(group)!=12 or len({r['network'] for r in group})!=12:raise ValueError('baseline group')
  ratios=[F(r['upper'])/((F(r['reference_lower'])+F(r['reference_upper']))/2) for r in group]
  comparisons.append({'method':method,'n':len(group),'median_upper_over_reference_midpoint':float(median(ratios)),'accept_1e-5':sum(r['accepted']['1/100000'] for r in group)})
 out={'extension':extension,'extension_records_sha256':sha(D/'records.json'),'baseline_64_cells':comparisons,'baseline_records_sha256':sha(BASE/'records.json'),'scope':'post-exit stored-value arithmetic; baseline 152 also independently checked full request manifest'}
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
