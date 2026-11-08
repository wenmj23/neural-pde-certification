"""Post-hoc 1%-margin tolerance analysis of saved certificate bounds."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,os,tempfile
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results';OUT=B/'2026-11-08/near_threshold_posthoc'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(l,u,t):return 'within' if u<=t else 'exceeds' if l>t else 'unknown'
def main():
 if OUT.exists():raise FileExistsError(OUT)
 rows=[];hashes={}
 for family in ('poly','sine'):
  p=B/f'2026-11-07/midterm_{family}_decisions/records.json';audit=json.loads((B/f'2026-11-07/midterm_{family}_decisions_readback.json').read_text())
  if sha(p)!=audit['records_sha256']:raise ValueError('source digest')
  hashes[family]=sha(p)
  for r in json.loads(p.read_text()):
   if r['status']!='complete':raise ValueError('incomplete source')
   l,u,rl,ru=map(F,(r['lower'],r['upper'],r['reference_lower'],r['reference_upper']));mid=(rl+ru)/2
   for factor in (F(99,100),F(101,100)):
    t=mid*factor;rows.append({'family':family,'network':r['network'],'degree':r['degree'],'mesh':r['mesh'],'factor':str(factor),'tolerance':str(t),'decision':dec(l,u,t),'reference_decision':dec(rl,ru,t)})
 OUT.parent.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix='.near_threshold_',dir=OUT.parent));(stage/'records.json').write_text(json.dumps(rows,indent=2)+'\n');s={'expected':432,'records':len(rows),'inputs_sha256':hashes,'records_sha256':sha(stage/'records.json'),'script_sha256':sha(Path(__file__)),'scope':'post-hoc tolerance = saved reference midpoint times 0.99 or 1.01; reclassification of existing ranges, no new certify calls; not blind, no timing'};(stage/'summary.json').write_text(json.dumps(s,indent=2)+'\n');os.rename(stage,OUT);print(json.dumps(s));return int(len(rows)!=432)
if __name__=='__main__':raise SystemExit(main())
