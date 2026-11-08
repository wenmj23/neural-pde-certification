"""Compare paired widths and describe an exact mesh-node aliasing family."""
from fractions import Fraction as F
from collections import Counter
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results/2026-11-07';SRC=B/'analytic_report_controls/records.json';OUT=B/'analytic_controls_pairwise_readback.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 if OUT.exists():raise FileExistsError(OUT)
 prior=json.loads((B/'analytic_report_controls_readback.json').read_text())
 if sha(SRC)!=prior['records_sha256']:raise ValueError('digest')
 rows=json.loads(SRC.read_text());idx={(r['family'],r['frequency'],r['mesh'],r['degree']):r for r in rows}
 if len(idx)!=78:raise ValueError('duplicate records')
 counts=Counter();quartic=[]
 for r in rows:
  if r['degree']!=2:continue
  key=(r['family'],r['frequency'],r['mesh'],3);other=idx[key];w2=F(r['upper'])-F(r['lower']);w3=F(other['upper'])-F(other['lower']);label='cubic_narrower' if w3<w2 else 'equal' if w3==w2 else 'cubic_wider';counts[(r['family'],label)]+=1
  if r['family']=='quartic':quartic.append({'mesh':r['mesh'],'quadratic_width':str(w2),'cubic_width':str(w3),'quadratic_width_approx':float(w2),'cubic_width_approx':float(w3)})
 if counts!={('sine','cubic_narrower'):36,('quartic','cubic_wider'):3}:raise ValueError('pairwise mismatch')
 # e=A sin(2*m*pi*x), nodes i/m: sin(2*i*pi)=0; at x=1/(4*m), sin(pi/2)=1.
 alias=[{'mesh':m,'frequency':2*m,'node_phase_factor':'2*i','unsampled_peak_x':str(F(1,4*m)),'exact_max_error':'1/10000'} for m in (8,16,32,64)]
 out={'paired_counts':[{'family':k[0],'comparison':k[1],'count':v} for k,v in sorted(counts.items())],'quartic':quartic,'alias_constructions':alias,'records_sha256':sha(SRC),'scope':'paired exact stored endpoints; alias cases follow analytic trigonometric identity and are manufactured, not trained models'}
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
