"""Post-hoc extension after 263's k=32 cases; not an independent holdout."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,os,tempfile
from arb_certification import exact_fraction
from arb_cubic_reconstruction import SineSix,CubicVerifier
from arb_shaped_remainder import ShapedVerifier
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results/2026-11-07';OUT=B/'analytic_sine32_extended';SRC=B/'analytic_report_controls'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 if OUT.exists():raise FileExistsError(OUT)
 prior=json.loads((SRC/'summary.json').read_text())
 if sha(SRC/'records.json')!=prior['records_sha256']:raise ValueError('prior digest')
 rows=[]
 for mesh in (512,1024):
  for degree in (2,3):
   p=SineSix(32,103);v=ShapedVerifier(p,'moment',subdivision=2) if degree==2 else CubicVerifier(p,'moment',subdivision=2);v.uniform(mesh);r=v.certify();L,U=map(exact_fraction,(r['lower'],r['upper']));K=F(1,10000)
   rows.append({'mesh':mesh,'degree':degree,'lower':str(L),'upper':str(U),'known_error':str(K),'contains':L<=K<=U,'ratio_met':L>0 and U<=F(51,50)*L,'ratio':str(U/L) if L>0 else None});print(mesh,degree,flush=True)
 stage=Path(tempfile.mkdtemp(prefix='.sine32_',dir=B));(stage/'records.json').write_text(json.dumps(rows,indent=2)+'\n');s={'records':len(rows),'contained':sum(x['contains'] for x in rows),'records_sha256':sha(stage/'records.json'),'prior_records_sha256':sha(SRC/'records.json'),'script_sha256':sha(Path(__file__)),'scope':'post-hoc analytic sine32 controls, no timing study'};(stage/'summary.json').write_text(json.dumps(s,indent=2)+'\n');os.rename(stage,OUT);print(json.dumps(s));return int(s['records']!=4 or s['contained']!=4)
if __name__=='__main__':raise SystemExit(main())
