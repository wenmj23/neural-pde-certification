"""Follow-up of high-frequency known-error sine tests with larger fixed grids."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json,os,sys
from arb_cubic_reconstruction import SineSix,CubicVerifier
from arb_shaped_remainder import ShapedVerifier
from arb_certification import exact_fraction
FREQ=(16,32); MESH=(128,256,512); DEGREE=(2,3)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 rows=[]
 for freq in FREQ:
  for n in MESH:
   for degree in DEGREE:
    p=SineSix(freq);v=ShapedVerifier(p,'moment') if degree==2 else CubicVerifier(p,'moment');v.uniform(n);r=v.certify();lo=exact_fraction(r['lower']);hi=exact_fraction(r['upper']);truth=F(1,10000)
    rows.append({'frequency':freq,'intervals':n,'degree':degree,'lower':str(lo),'upper':str(hi),'known_exact_error':str(truth),'contains_truth':lo<=truth<=hi,'gap_goal_1_02':lo>0 and hi<=F(51,50)*lo})
    print(freq,n,degree,'contains',lo<=truth<=hi,flush=True)
 data=(json.dumps(rows,indent=2)+'\n').encode();tmp=a.out/'records.tmp'
 with tmp.open('xb') as f:f.write(data);f.flush();os.fsync(f.fileno())
 os.replace(tmp,a.out/'records.json');s={'frequencies':FREQ,'meshes':MESH,'degrees':DEGREE,'expected':len(FREQ)*len(MESH)*len(DEGREE),'records':len(rows),'contains_truth':sum(r['contains_truth'] for r in rows),'goal_met':sum(r['gap_goal_1_02'] for r in rows),'sha256':hashlib.sha256(data).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'failures':[r for r in rows if not r['contains_truth']],'scope':'synthetic analytic sine errors, not trained PINNs; conditional on interval implementation'}
 (a.out/'summary.json').write_text(json.dumps(s,indent=2)+'\n');print(json.dumps(s));return int(bool(s['failures']))
if __name__=='__main__':sys.exit(main())
