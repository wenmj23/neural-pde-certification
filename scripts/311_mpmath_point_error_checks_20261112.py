"""Interval enclosures of error at exact rational points with mpmath.iv."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,math,os,tempfile
import mpmath
from interval_certification import context,endpoint_fraction
ROOT=Path(__file__).resolve().parents[1]
CFG=ROOT/'configs/mpmath_point_error_checks_2026-11-12.json'
DEST=ROOT/'results/2026-11-12/mpmath_point_error_checks'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def interval_candidate(net,iv,x):
    def exact(v):
        if not isinstance(v,(float,int)) or isinstance(v,bool) or not math.isfinite(v):raise ValueError('numeric input')
        z=iv.mpf(float(v)); q=F.from_float(float(v))
        if not endpoint_fraction(z)==endpoint_fraction(z,1)==q:raise ArithmeticError('weight conversion')
        return z
    a=[x]
    for layer in net['layers']:
        w=[[exact(q) for q in row] for row in layer['weight']];bias=[exact(q) for q in layer['bias']]
        z=[sum((v*t for v,t in zip(row,a)),b) for row,b in zip(w,bias)]
        if layer['activation']=='tanh':a=[1-2/(iv.exp(2*t)+1) for t in z]
        elif layer['activation']=='identity':a=z
        else:raise ValueError('activation')
    if len(a)!=1:raise ValueError('output')
    u=sum((iv.mpf(c)*(x-x**(j+2))/((j+1)*(j+2)) for j,c in enumerate(net['source_coefficients'])),iv.mpf(0))
    return x*(1-x)*a[0]-u

def main():
    if DEST.exists():raise FileExistsError(DEST)
    c=json.loads(CFG.read_text());inputs={k:ROOT/c[k+'_path'] for k in ['fixture','baseline','spots']}
    for k,p in inputs.items():
        if sha(p)!=c[k+'_sha256']:raise ValueError('input hash')
    nets={n['label']:n for n in json.loads(inputs['fixture'].read_text())['networks']}
    old={(r['network'],r['method']):r for r in json.loads(inputs['baseline'].read_text()) if r['intervals']==c['certificate_mesh']}
    spots=json.loads(inputs['spots'].read_text())
    if len(nets)!=12 or len(old)!=48 or len(spots)!=48:raise ValueError('manifest')
    rows=[]
    for spot in spots:
        name=spot['network'];iv=context(c['digits']);x=iv.mpf(spot['grid_index'])/spot['grid_intervals']
        err=interval_candidate(nets[name],iv,x);lo,hi=endpoint_fraction(err),endpoint_fraction(err,1)
        lower=max(F(0),lo if lo>0 else (-hi if hi<0 else F(0)));upper=max(abs(lo),abs(hi))
        point_num=F(spot['error_110']);bounds={m:F(old[(name,m)]['upper']) for m in ['signed_quadratic','signed_cubic']}
        rows.append({'network':name,'kind':spot['kind'],'point':str(F(spot['grid_index'],spot['grid_intervals'])),
          'signed_error_lower':str(lo),'signed_error_upper':str(hi),'absolute_error_lower':str(lower),'absolute_error_upper':str(upper),
          'high_precision_error':str(point_num),'contains_high_precision_error':lo<=point_num<=hi,
          'saved_upper':{k:str(v) for k,v in bounds.items()},
          'point_proved_below_saved_upper':{k:upper<=v for k,v in bounds.items()},
          'point_proved_above_saved_upper':{k:lower>v for k,v in bounds.items()}})
    if len(rows)!=48:raise ValueError('output count')
    DEST.parent.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix='.point_error_iv_',dir=DEST.parent))
    (stage/'locked_config.json').write_bytes(CFG.read_bytes());(stage/'records.json').write_text(json.dumps(rows,indent=2)+'\n')
    out={'expected':48,'records':len(rows),'inputs':{k:sha(p) for k,p in inputs.items()},'config_sha256':sha(CFG),
         'script_sha256':sha(Path(__file__)),'helper_sha256':sha(ROOT/'scripts/interval_certification.py'),
         'records_sha256':sha(stage/'records.json'),'mpmath_version':mpmath.__version__,'scope':c['scope']}
    (stage/'summary.json').write_text(json.dumps(out,indent=2)+'\n');os.rename(stage,DEST);print(json.dumps(out))
if __name__=='__main__':main()
