"""Numerical root and Green checks near previously selected grid peaks."""
from fractions import Fraction
from pathlib import Path
import hashlib,importlib.util,json,os,tempfile,time
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1]
CFG=ROOT/'configs/polynomial_peak_audit_2026-11-12.json'
DEST=ROOT/'results/2026-11-12/polynomial_peak_audit'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
spec=importlib.util.spec_from_file_location('quad_model',ROOT/'scripts/296_green_quadrature_audit_20261112.py')
helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)

def peak(net,spot,digits):
    with mp.workdps(digits):
        r,e=helper.model(net);h=mp.mpf(1)/spot['grid_intervals'];guess=spot['grid_index']*h
        x=mp.findroot(lambda z:mp.diff(e,z),(guess-h,guess+h),tol=mp.mpf(10)**(-digits+15),maxsteps=50)
        value=e(x);slope=mp.diff(e,x);curvature=-r(x)
        if not 0<x<1:raise ValueError('root outside domain')
        return {'digits':digits,'root':mp.nstr(x,digits-6),'error':mp.nstr(value,digits-6),
                'abs_derivative':mp.nstr(abs(slope),8),'curvature':mp.nstr(curvature,digits-6),
                'local_abs_max_curvature_test':bool(value*curvature<0)}

def main():
    if DEST.exists():raise FileExistsError(DEST)
    c=json.loads(CFG.read_text()); paths={k:ROOT/c[k+'_path'] for k in ['fixture','references','spots']}
    for k,p in paths.items():
        if sha(p)!=c[k+'_sha256']:raise ValueError('input hash '+k)
    nets=json.loads(paths['fixture'].read_text())['networks'];refs={r['network']:r for r in json.loads(paths['references'].read_text())}
    spots={r['network']:r for r in json.loads(paths['spots'].read_text()) if r['kind']=='max_abs_error'}
    if len(nets)!=12 or set(refs)!=set(spots) or set(refs)!={n['label'] for n in nets}:raise ValueError('manifest')
    DEST.parent.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix='.peak_audit_',dir=DEST.parent))
    (stage/'locked_config.json').write_bytes(CFG.read_bytes());rows=[];start=time.monotonic()
    for i,n in enumerate(nets):
        label=n['label'];roots=[peak(n,spots[label],d) for d in c['root_digits']]
        p=roots[-1]['root'];quad=[helper.evaluate(n,p,d) for d in c['quadrature_digits']]
        ref=refs[label];lo,hi=Fraction(ref['lower']),Fraction(ref['upper'])
        a=abs(Fraction(roots[-1]['error']))
        row={'network':label,'initial_grid_index':spots[label]['grid_index'],
             'initial_grid_intervals':spots[label]['grid_intervals'], 'roots':roots,'quadrature':quad,
             'reference_lower':str(lo),'reference_upper':str(hi),'reference_converged':ref['converged'],
             'numerical_peak_within_reference':lo<=a<=hi}
        (stage/f'case_{i:02d}.json').write_text(json.dumps(row,indent=2)+'\n');rows.append(row)
        print(label,'complete',flush=True)
    if len(rows)!=12:raise ValueError('result count')
    (stage/'records.json').write_text(json.dumps(rows,indent=2)+'\n')
    summary={'expected':12,'records':len(rows),'inputs':{k:sha(p) for k,p in paths.items()},
             'script_sha256':sha(Path(__file__)),'helper_sha256':sha(Path(helper.__file__)),
             'config_sha256':sha(CFG),'records_sha256':sha(stage/'records.json'),
             'cases':{p.name:sha(p) for p in sorted(stage.glob('case_*.json'))},
             'elapsed_seconds':time.monotonic()-start,'scope':c['scope']}
    (stage/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');os.rename(stage,DEST)
    print(json.dumps({k:v for k,v in summary.items() if k!='cases'}),flush=True)
if __name__=='__main__':main()
