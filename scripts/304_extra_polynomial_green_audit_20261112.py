"""Green and source-omission controls on 24 additional old polynomial PINNs."""
from fractions import Fraction
from pathlib import Path
import hashlib,importlib.util,json,os,tempfile,time
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1]
CFG=ROOT/'configs/extra_polynomial_green_audit_2026-11-12.json'
DEST=ROOT/'results/2026-11-12/extra_polynomial_green_audit'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
sp=importlib.util.spec_from_file_location('quad_model',ROOT/'scripts/296_green_quadrature_audit_20261112.py')
helper=importlib.util.module_from_spec(sp);sp.loader.exec_module(helper)

def source_control(net,point):
    with mp.workdps(80):
        r,e=helper.model(net);q=Fraction(point);x=mp.mpf(q.numerator)/q.denominator
        direct_r=r(x);independent_r=-mp.diff(e,x,2)
        coeff=[mp.mpf(z) for z in net['source_coefficients']]
        f=sum((c*x**j for j,c in enumerate(coeff)),mp.mpf(0))
        injected=r(x)+f
        return {'correct_residual':mp.nstr(direct_r,72),
                'independent_minus_error_second_derivative':mp.nstr(independent_r,72),
                'absolute_disagreement':mp.nstr(abs(direct_r-independent_r),10),
                'source_value':mp.nstr(f,72),
                'injected_omission_residual':mp.nstr(injected,72),
                'injected_fault_nonzero':bool(abs(injected-independent_r)>mp.mpf('1e-50'))}

def main():
    if DEST.exists():raise FileExistsError(DEST)
    c=json.loads(CFG.read_text());groups=[];labels=[]
    for g in c['groups']:
        fixture,old=ROOT/g['fixture_path'],ROOT/g['records_path']
        if sha(fixture)!=g['fixture_sha256'] or sha(old)!=g['records_sha256']:raise ValueError('input hash')
        nets=json.loads(fixture.read_text())['networks'];prior=json.loads(old.read_text())
        idx={r['network']:r for r in prior if r['degree']==c['degree'] and r['intervals']==c['mesh']}
        if len(nets)!=12 or len(idx)!=12 or {n['label'] for n in nets}!=set(idx):raise ValueError('manifest')
        groups.append((g,nets,idx));labels.extend(n['label'] for n in nets)
    if len(set(labels))!=24:raise ValueError('network names')
    DEST.parent.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix='.extra_poly_green_',dir=DEST.parent))
    (stage/'locked_config.json').write_bytes(CFG.read_bytes());rows=[];t0=time.monotonic();ni=0
    for g,nets,idx in groups:
        for n in nets:
            upper=Fraction(idx[n['label']]['upper'])
            for pi,point in enumerate(c['points']):
                vals=[helper.evaluate(n,point,d) for d in c['digits']]
                control=source_control(n,point)
                with mp.workdps(100):
                    u=mp.mpf(upper.numerator)/upper.denominator
                    within=bool(abs(mp.mpf(vals[-1]['direct_error']))<=u)
                row={'group':g['name'],'network':n['label'],'source':n['source_name'],
                     'point':point,'results':vals,'source_control':control,
                     'saved_cubic_upper':str(upper),'direct_below_saved_cubic_upper':within}
                (stage/f'case_{ni:02d}_{pi:02d}.json').write_text(json.dumps(row,indent=2)+'\n');rows.append(row)
                print(n['label'],point,'complete',flush=True)
            ni+=1
    if len(rows)!=72:raise ValueError('result count')
    (stage/'records.json').write_text(json.dumps(rows,indent=2)+'\n')
    summary={'expected':72,'records':len(rows),'config_sha256':sha(CFG),
       'script_sha256':sha(Path(__file__)),'helper_sha256':sha(Path(helper.__file__)),
       'records_sha256':sha(stage/'records.json'),
       'case_sha256':{p.name:sha(p) for p in sorted(stage.glob('case_*.json'))},
       'elapsed_seconds':time.monotonic()-t0,'scope':c['scope']}
    if len(summary['case_sha256'])!=72:raise ValueError('case count')
    (stage/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');os.rename(stage,DEST)
    print('COMPLETE',json.dumps({k:v for k,v in summary.items() if k!='case_sha256'}),flush=True)
if __name__=='__main__':main()
