"""Alternative mpmath.iv global residual bounds for frozen polynomial PINNs.

For each cell, bound the signed residual directly by interval forward-mode
derivatives and take max |r|/8. No neural training and no Arb calls.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,math,os,tempfile,time
import mpmath
from interval_certification import context,endpoint_fraction
ROOT=Path(__file__).resolve().parents[1]
CFG=ROOT/'configs/mpmath_global_residual_audit_2026-11-12.json'
DEST=ROOT/'results/2026-11-12/mpmath_global_residual_audit'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def interval_model(net,iv):
    conversions=0
    def exact(value):
        nonlocal conversions
        if not isinstance(value,(int,float)) or isinstance(value,bool) or not math.isfinite(value):raise ValueError('bad numeric input')
        v=float(value);z=iv.mpf(v);q=F.from_float(v)
        if not endpoint_fraction(z)==endpoint_fraction(z,1)==q:raise ArithmeticError('non-exact weight import')
        conversions+=1;return z
    layers=[([[exact(v) for v in row] for row in layer['weight']],
             [exact(v) for v in layer['bias']],layer['activation']) for layer in net['layers']]
    # These source coefficients are small integers in this experiment.
    if any(not isinstance(v,int) or isinstance(v,bool) for v in net['source_coefficients']):raise ValueError('noninteger source')
    coeff=[iv.mpf(v) for v in net['source_coefficients']]
    def residual(x):
        a=[x];a1=[iv.mpf(1)];a2=[iv.mpf(0)]
        for w,b,act in layers:
            z=[sum((v*q for v,q in zip(row,a)),bias) for row,bias in zip(w,b)]
            z1=[sum((v*q for v,q in zip(row,a1)),iv.mpf(0)) for row in w]
            z2=[sum((v*q for v,q in zip(row,a2)),iv.mpf(0)) for row in w]
            if act=='tanh':
                a=[1-2/(iv.exp(2*q)+1) for q in z]
                s=[1-q**2 for q in a]
                a1=[v*u for v,u in zip(s,z1)]
                a2=[v*u-2*t*v*d**2 for v,u,t,d in zip(s,z2,a,z1)]
            elif act=='identity':a,a1,a2=z,z1,z2
            else:raise ValueError('unknown activation')
        if len(a)!=1:raise ValueError('network output')
        f=sum((c*x**j for j,c in enumerate(coeff)),iv.mpf(0))
        return -x*(1-x)*a2[0]-2*(1-2*x)*a1[0]+2*a[0]-f
    return residual,conversions

def main():
    if DEST.exists():raise FileExistsError(DEST)
    c=json.loads(CFG.read_text());inputs={k:ROOT/c[k+'_path'] for k in ['fixture','baseline','spots']}
    for k,p in inputs.items():
        if sha(p)!=c[k+'_sha256']:raise ValueError('input hash '+k)
    nets=json.loads(inputs['fixture'].read_text())['networks'];base=json.loads(inputs['baseline'].read_text());sp=json.loads(inputs['spots'].read_text())
    old={(r['network'],r['method']):r for r in base if r['intervals']==64}
    if len(nets)!=12 or len(old)!=48 or len(sp)!=48:raise ValueError('input manifest')
    DEST.parent.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix='.iv_global_',dir=DEST.parent))
    (stage/'locked_config.json').write_bytes(CFG.read_bytes());rows=[];pointchecks=[];t0=time.monotonic()
    for ni,net in enumerate(nets):
        iv=context(c['digits']);residual,conversions=interval_model(net,iv);name=net['label']
        for spot in [s for s in sp if s['network']==name]:
            x=iv.mpf(spot['grid_index'])/spot['grid_intervals'];v=residual(x)
            lo,hi=endpoint_fraction(v),endpoint_fraction(v,1);hp=F(spot['residual_110'])
            pointchecks.append({'network':name,'kind':spot['kind'],'lower':str(lo),'upper':str(hi),
                                'high_precision_residual':str(hp),'contains_high_precision_point':lo<=hp<=hi})
        ref=old[(name,'global')];rl,ru=F(ref['reference_lower']),F(ref['reference_upper']);old_global=F(ref['upper'])
        for mesh in c['meshes']:
            max_r=F(0);max_cell=None;start=time.monotonic()
            for i in range(mesh):
                cell=iv.mpf([iv.mpf(i)/mesh,iv.mpf(i+1)/mesh]);r=residual(cell)
                lo,hi=endpoint_fraction(r),endpoint_fraction(r,1);candidate=max(abs(lo),abs(hi))
                if candidate>max_r:max_r,max_cell=candidate,i
            U=max_r/8
            row={'network':name,'mesh':mesh,'precision_digits':c['digits'],'exact_weight_imports':conversions,
                 'max_residual_upper':str(max_r),'upper':str(U),'max_cell':max_cell,
                 'old_arb_global_upper':str(old_global),'reference_lower':str(rl),'reference_upper':str(ru),
                 'contains_saved_reference':ru<=U,'upper_over_reference_midpoint':float(U/((rl+ru)/2)),
                 'elapsed_seconds':time.monotonic()-start}
            (stage/f'case_{ni:02d}_{mesh}.json').write_text(json.dumps(row,indent=2)+'\n');rows.append(row)
            print(name,mesh,'cells complete, seconds',round(row['elapsed_seconds'],2),flush=True)
    if len(rows)!=36 or len(pointchecks)!=48:raise ValueError('result count')
    (stage/'records.json').write_text(json.dumps(rows,indent=2)+'\n')
    (stage/'point_checks.json').write_text(json.dumps(pointchecks,indent=2)+'\n')
    summary={'expected':36,'records':len(rows),'expected_point_checks':48,'point_checks':len(pointchecks),
      'inputs':{k:sha(p) for k,p in inputs.items()},'config_sha256':sha(CFG),
      'script_sha256':sha(Path(__file__)),'interval_helper_sha256':sha(ROOT/'scripts/interval_certification.py'),
      'records_sha256':sha(stage/'records.json'),'points_sha256':sha(stage/'point_checks.json'),
      'cases':{p.name:sha(p) for p in sorted(stage.glob('case_*.json'))},
      'elapsed_seconds':time.monotonic()-t0,'mpmath_version':mpmath.__version__,'scope':c['scope']}
    if len(summary['cases'])!=36:raise ValueError('case count')
    (stage/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');os.rename(stage,DEST)
    print('COMPLETE',json.dumps({k:v for k,v in summary.items() if k!='cases'}),flush=True)
if __name__=='__main__':main()
