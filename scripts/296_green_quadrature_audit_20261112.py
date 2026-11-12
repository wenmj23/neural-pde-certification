"""Independent high-precision pointwise Green-integral diagnostics.

Residuals use a separate manual forward derivative implementation; exact
polynomial solutions use the known forcing. mpmath.quad is not interval proof.
Each case is staged separately, so a failed or interrupted run is visible.
"""
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import os
import platform
import sys
import tempfile
import time

import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
CFG=ROOT/'configs/green_quadrature_audit_2026-11-12.json'
DEST=ROOT/'results/2026-11-12/green_quadrature_audit'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def model(net):
    layers=[(tuple(tuple(mp.mpf(v) for v in row) for row in layer['weight']),
             tuple(mp.mpf(v) for v in layer['bias']),layer['activation'])
            for layer in net['layers']]
    coeff=tuple(mp.mpf(v) for v in net['source_coefficients'])
    def candidate(x, derivatives=False):
        a=[x];da=[mp.mpf(1)];dda=[mp.mpf(0)]
        for weights,bias,activation in layers:
            z=[sum((w*q for w,q in zip(row,a)),b) for row,b in zip(weights,bias)]
            if derivatives:
                dz=[sum((w*q for w,q in zip(row,da)),mp.mpf(0)) for row in weights]
                ddz=[sum((w*q for w,q in zip(row,dda)),mp.mpf(0)) for row in weights]
            if activation=='tanh':
                a=[mp.tanh(q) for q in z]
                if derivatives:
                    da=[(1-v*v)*v1 for v,v1 in zip(a,dz)]
                    dda=[(1-v*v)*v2-2*v*(1-v*v)*v1*v1 for v,v1,v2 in zip(a,dz,ddz)]
            elif activation=='identity':
                a=z
                if derivatives:da,dda=dz,ddz
            else:raise ValueError('activation')
        if len(a)!=1:raise ValueError('output')
        if derivatives:return x*(1-x)*dda[0]+2*(1-2*x)*da[0]-2*a[0]
        return x*(1-x)*a[0]
    def rhs(x):return sum((c*x**j for j,c in enumerate(coeff)),mp.mpf(0))
    def exact(x):return sum((c*(x-x**(j+2))/((j+1)*(j+2)) for j,c in enumerate(coeff)),mp.mpf(0))
    def residual(x):return -candidate(x,True)-rhs(x)
    def error(x):return candidate(x)-exact(x)
    return residual,error

def evaluate(net, point, digits):
    start=time.monotonic()
    with mp.workdps(digits):
        residual,error=model(net)
        f=Fraction(point);x=mp.mpf(f.numerator)/f.denominator
        a,ae=mp.quad(lambda s:s*residual(s),[0,x],error=True)
        b,be=mp.quad(lambda s:(1-s)*residual(s),[x,1],error=True)
        green=(1-x)*a+x*b
        direct=error(x)
        return {'digits':digits,'green':mp.nstr(green,digits-6),
                'direct_error':mp.nstr(direct,digits-6),
                'difference':mp.nstr(abs(green-direct),10),
                'quad_error_estimate':mp.nstr((1-x)*ae+x*be,10),
                'elapsed_seconds':time.monotonic()-start}

def main():
    if DEST.exists():raise FileExistsError(DEST)
    c=json.loads(CFG.read_text());fixture=ROOT/c['fixture_path'];baseline=ROOT/c['baseline_path']
    if sha(fixture)!=c['fixture_sha256'] or sha(baseline)!=c['baseline_sha256']:raise ValueError('input hash')
    nets=json.loads(fixture.read_text())['networks'];old=json.loads(baseline.read_text())
    idx={(r['network'],r['method']):r for r in old if r['intervals']==c['grid_intervals_for_saved_certificate']}
    if len(nets)!=12 or len({n['label'] for n in nets})!=12:raise ValueError('manifest')
    if any((n['label'],m) not in idx for n in nets for m in ('global','spatial','signed_quadratic','signed_cubic')):raise ValueError('baseline manifest')
    DEST.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix='.green_quad_stage_',dir=DEST.parent))
    (stage/'locked_config.json').write_bytes(CFG.read_bytes())
    (stage/'in_progress.json').write_text(json.dumps({'expected':len(nets)*len(c['points']),'scope':c['scope']},indent=2)+'\n')
    rows=[];start=time.monotonic()
    for ni,net in enumerate(nets):
        name=net['label'];upper=Fraction(idx[(name,'signed_cubic')]['upper'])
        for pi,point in enumerate(c['points']):
            values=[evaluate(net,point,d) for d in c['digits']]
            with mp.workdps(100):
                discrepancy=abs(mp.mpf(values[0]['green'])-mp.mpf(values[1]['green']))
                direct=abs(mp.mpf(values[1]['direct_error']))
                up=mp.mpf(upper.numerator)/upper.denominator
                row={'network':name,'point':point,'results':values,
                     'green_cross_precision_difference':mp.nstr(discrepancy,12),
                     'direct_below_saved_cubic_upper':bool(direct<=up),
                     'saved_cubic_upper_fraction':str(upper)}
            p=stage/f'case_{ni:02d}_{pi:02d}.json'
            p.write_text(json.dumps(row,indent=2,allow_nan=False)+'\n')
            rows.append(row)
            print(name,point,'complete; seconds',round(sum(v['elapsed_seconds'] for v in values),2),flush=True)
    (stage/'records.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
    s={'expected':84,'records':len(rows),'fixture_sha256':sha(fixture),'baseline_sha256':sha(baseline),
       'script_sha256':sha(Path(__file__)),'config_sha256':sha(CFG),'records_sha256':sha(stage/'records.json'),
       'case_sha256':{p.name:sha(p) for p in sorted(stage.glob('case_*.json'))},
       'elapsed_seconds':time.monotonic()-start,'mpmath_version':mp.__version__,
       'python':sys.version,'platform':platform.platform(),'scope':c['scope']}
    if len(rows)!=84 or len(s['case_sha256'])!=84:raise ValueError('result count')
    (stage/'summary.json').write_text(json.dumps(s,indent=2,allow_nan=False)+'\n')
    os.rename(stage,DEST)
    print('ALL CASES COMPLETE',json.dumps({k:v for k,v in s.items() if k!='case_sha256'}),flush=True)

if __name__=='__main__':main()
