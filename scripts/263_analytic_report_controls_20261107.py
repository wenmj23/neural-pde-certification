"""Manufactured controls with analytically known continuous maximum error."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib,json,os,tempfile
from flint import arb
from arb_certification import exact_fraction
from arb_problem_families import QuarticError
from arb_cubic_reconstruction import SineSix,CubicVerifier
from arb_shaped_remainder import ShapedVerifier

ROOT=Path(__file__).resolve().parents[1]
CONFIG=ROOT/'configs/analytic_report_controls_2026-11-07.json'
OUT=ROOT/'results/2026-11-07/analytic_report_controls'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 if OUT.exists():raise FileExistsError(OUT)
 c=json.loads(CONFIG.read_text()); known=F(c['absolute_exact_error']);records=[]
 class QuarticSix(QuarticError):
  def r_fourth(self,x):self.guard();return arb(0)
 requests=[('sine',k,m,d) for k,m,d in product(c['frequencies'],c['meshes'],c['degrees'])]+[('quartic',None,m,d) for m,d in product(c['quartic_meshes'],c['degrees'])]
 for family,k,mesh,degree in requests:
  problem=SineSix(k,c['precision_bits']) if family=='sine' else QuarticSix(c['precision_bits'])
  v=ShapedVerifier(problem,'moment',subdivision=c['bernstein_subdivision_depth']) if degree==2 else CubicVerifier(problem,'moment',subdivision=c['bernstein_subdivision_depth'])
  v.uniform(mesh);z=v.certify();L,U=map(exact_fraction,(z['lower'],z['upper']));D=exact_fraction(z['remainder_upper'])
  decisions={t:('within' if U<=F(t) else 'exceeds' if L>F(t) else 'undetermined') for t in c['tolerances']}
  records.append({'family':family,'frequency':k,'mesh':mesh,'degree':degree,'lower':str(L),'upper':str(U),'remainder_upper':str(D),'known_error':str(known),'contains_exact':L<=known<=U,'ratio_met':L>0 and U<=F(c['ratio_target'])*L,'decisions':decisions})
  print(family,k,mesh,degree,flush=True)
 OUT.parent.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix='.analytic_report_',dir=OUT.parent));(stage/'locked_config.json').write_bytes(CONFIG.read_bytes());(stage/'records.json').write_text(json.dumps(records,indent=2)+'\n')
 s={'expected':len(requests),'records':len(records),'contains_exact':sum(r['contains_exact'] for r in records),'records_sha256':sha(stage/'records.json'),'config_sha256':sha(CONFIG),'script_sha256':sha(Path(__file__)),'scope':c['scope']};(stage/'summary.json').write_text(json.dumps(s,indent=2)+'\n');os.rename(stage,OUT);print(json.dumps(s));return int(s['records']!=78 or s['contains_exact']!=78)
if __name__=='__main__':raise SystemExit(main())
