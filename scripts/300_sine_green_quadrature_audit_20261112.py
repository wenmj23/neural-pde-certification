"""Finite-point sine-forcing Green identity check, not interval certification."""
from fractions import Fraction
from pathlib import Path
import hashlib,importlib.util,json,os,platform,sys,tempfile,time
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1]
CFG=ROOT/'configs/sine_green_quadrature_audit_2026-11-12.json'
DEST=ROOT/'results/2026-11-12/sine_green_quadrature_audit'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
spec=importlib.util.spec_from_file_location('quad_model',ROOT/'scripts/296_green_quadrature_audit_20261112.py')
helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)

def evaluate(n,point,digits):
    start=time.monotonic()
    with mp.workdps(digits):
        z=dict(n);z['source_coefficients']=[]
        neg_v2,v=helper.model(z)
        k=mp.mpf(n['frequency']);freq=k*mp.pi
        def r(x):return neg_v2(x)-freq**2*mp.sin(freq*x)
        def e(x):return v(x)-mp.sin(freq*x)
        a=Fraction(point);x=mp.mpf(a.numerator)/a.denominator
        left,le=mp.quad(lambda s:s*r(s),[0,x],error=True)
        right,re=mp.quad(lambda s:(1-s)*r(s),[x,1],error=True)
        green=(1-x)*left+x*right;direct=e(x)
        return {'digits':digits,'green':mp.nstr(green,digits-6),
                'direct_error':mp.nstr(direct,digits-6),
                'difference':mp.nstr(abs(green-direct),10),
                'quad_error_estimate':mp.nstr((1-x)*le+x*re,10),
                'elapsed_seconds':time.monotonic()-start}

def main():
    if DEST.exists():raise FileExistsError(DEST)
    c=json.loads(CFG.read_text());fixture=ROOT/c['fixture_path'];oldpath=ROOT/c['records_path']
    if sha(fixture)!=c['fixture_sha256'] or sha(oldpath)!=c['records_sha256']:raise ValueError('input hash')
    nets=json.loads(fixture.read_text())['networks'];prior=json.loads(oldpath.read_text())
    idx={r['network']:r for r in prior if r['degree']==c['certificate_degree'] and r['intervals']==c['certificate_mesh']}
    if len(nets)!=12 or len(idx)!=12 or {n['label'] for n in nets}!=set(idx):raise ValueError('manifest')
    DEST.parent.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix='.sine_green_stage_',dir=DEST.parent))
    (stage/'locked_config.json').write_bytes(CFG.read_bytes());rows=[];start=time.monotonic()
    for ni,n in enumerate(nets):
        upper=Fraction(idx[n['label']]['upper'])
        for pi,point in enumerate(c['points']):
            vals=[evaluate(n,point,d) for d in c['digits']]
            with mp.workdps(100):
                diff=abs(mp.mpf(vals[0]['green'])-mp.mpf(vals[1]['green']))
                u=mp.mpf(upper.numerator)/upper.denominator
                row={'network':n['label'],'frequency':n['frequency'],'point':point,'results':vals,
                     'green_cross_precision_difference':mp.nstr(diff,12),
                     'direct_below_saved_cubic_upper':bool(abs(mp.mpf(vals[-1]['direct_error']))<=u),
                     'saved_cubic_upper_fraction':str(upper)}
            (stage/f'case_{ni:02d}_{pi:02d}.json').write_text(json.dumps(row,indent=2)+'\n')
            rows.append(row);print(n['label'],point,'seconds',round(sum(v['elapsed_seconds'] for v in vals),2),flush=True)
    if len(rows)!=60:raise ValueError('record count')
    (stage/'records.json').write_text(json.dumps(rows,indent=2)+'\n')
    s={'expected':60,'records':len(rows),'fixture_sha256':sha(fixture),'old_records_sha256':sha(oldpath),
       'script_sha256':sha(Path(__file__)),'helper_sha256':sha(Path(helper.__file__)),
       'config_sha256':sha(CFG),'records_sha256':sha(stage/'records.json'),
       'case_sha256':{p.name:sha(p) for p in sorted(stage.glob('case_*.json'))},
       'elapsed_seconds':time.monotonic()-start,'python':sys.version,'platform':platform.platform(),
       'mpmath_version':mp.__version__,'scope':c['scope']}
    if len(s['case_sha256'])!=60:raise ValueError('case count')
    (stage/'summary.json').write_text(json.dumps(s,indent=2)+'\n');os.rename(stage,DEST)
    print('COMPLETE',json.dumps({k:v for k,v in s.items() if k!='case_sha256'}),flush=True)
if __name__=='__main__':main()
