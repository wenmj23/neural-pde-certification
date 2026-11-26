"""Locked fixed-mesh bound comparison and actual online signed-certifier runs.

No reference/true solution is read here. Evaluate serially for interpretable
within-process timing; retain all traces, unresolved cases and budget exits.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,importlib.util,json,os,tempfile,time,platform
import flint
from arb_certification import exact_fraction
from arb_cubic_reconstruction import SixthOrderDeep,CubicVerifier
from arb_shaped_remainder import ShapedVerifier
from arb_polynomial_forcing import PolynomialForcing

ROOT=Path(__file__).resolve().parents[1]
def module(name,path):
    sp=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
baseline=module('baseline151','scripts/151_new_pinn_matched_baselines_20261102.py').baseline
SplitProblem=module('split225','scripts/225_derivative_split_20261105.py').SplitProblem
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(path,obj):
    data=(json.dumps(obj,indent=2,allow_nan=False)+'\n').encode()
    with path.open('xb') as f:f.write(data);f.flush();os.fsync(f.fileno())
def construct(fixture,h,net,degree,pieces,c):
    if net['source_kind']=='sine':p=SixthOrderDeep(fixture,net['label'],h,c['precision_bits'])
    else:p=PolynomialForcing(fixture,net['label'],h,net['source_coefficients'],c['precision_bits'])
    p=SplitProblem(p,pieces)
    v=ShapedVerifier(p,'moment',subdivision=c['bernstein_depth']) if degree==2 else CubicVerifier(p,'moment',subdivision=c['bernstein_depth'])
    return p,v
def bounds(r):return {k:str(exact_fraction(r[k])) for k in ['lower','upper','z_lower','z_upper','remainder_upper']}
def decision(b,g):
    lo,hi=F(b['lower']),F(b['upper'])
    if g['goal']=='enclosure':
        met=lo==hi==0 or (lo>0 and hi<=F(g['gap_ratio'])*lo)
        return 'gap_met' if met else 'gap_unresolved',bool(met)
    tol=F(g['tolerance']);d='within' if hi<=tol else 'exceeds' if lo>tol else 'unknown'
    return d,d!='unknown'
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--training',type=Path,required=True);ap.add_argument('--config',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    c=json.loads(a.config.read_text());fixture=(a.training/'frozen_networks.json').resolve();h=sha(fixture);ts=json.loads((a.training/'summary.json').read_text())
    if h!=ts['file_sha256']['frozen_networks.json']:raise ValueError('training fixture hash')
    nets=json.loads(fixture.read_text())['networks'];a.out.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix='.fixed-evaluation-',dir=a.out.parent));write(stage/'locked_config.json',c)
    fixed=[];online=[]
    for ni,net in enumerate(nets):
        for degree in [2,3]:
            p,v=construct(fixture,h,net,degree,1,c);v.uniform(c['fixed_mesh']);r=v.certify();b=bounds(r)
            fixed.append({'network':net['label'],'method':'signed_quadratic' if degree==2 else 'signed_cubic','mesh':c['fixed_mesh'],**b})
            if degree==2:
                glob,spatial=baseline(v.prepare(v.nodes),v.depth)
                for method,u in [('global',glob),('spatial',spatial)]:fixed.append({'network':net['label'],'method':method,'mesh':c['fixed_mesh'],'lower':'0','upper':str(u)})
        write(stage/f'fixed_{ni:02d}.json',fixed[-4:])
        for si,strategy in enumerate(c['strategies']):
            for gi,goal in enumerate(c['goals']):
                for repeat in range(c['timing_repeats']):
                    tick=time.perf_counter();trace=[];last=None;state='unknown';met=False;stop='mesh_limit'
                    p,v=construct(fixture,h,net,strategy['degree'],strategy['derivative_splits'],c)
                    for mesh in c['meshes']:
                        if time.perf_counter()-tick>=c['soft_seconds_per_request']:stop='time_limit';break
                        start=time.perf_counter()
                        try:
                            v.uniform(mesh);r=v.certify();last=bounds(r);state,met=decision(last,goal)
                            attempt={'mesh':mesh,'status':'finite',**last,'decision':state,'goal_met':met}
                        except ArithmeticError as exc:
                            attempt={'mesh':mesh,'status':'arithmetic_unresolved','error':str(exc)}
                        attempt['seconds']=time.perf_counter()-start;trace.append(attempt)
                        if met:stop='goal_met';break
                        if time.perf_counter()-tick>=c['soft_seconds_per_request']:stop='time_limit';break
                    elapsed=time.perf_counter()-tick
                    row={'network':net['label'],'strategy':strategy,'goal':goal,'repeat':repeat,'stop_reason':stop,'goal_met':met,'decision':state,'enclosure':last,'trace':trace,'request_seconds':elapsed,'residual_evaluations':v.residual_evaluations,'derivative_interval_requests':v.curvature_evaluations,'underlying_derivative_calls':p.underlying_derivative_calls,'timing_scope':'constructor including frozen-file I/O and all attempts; excludes JSON export and outer process startup'}
                    online.append(row);write(stage/f'online_{ni:02d}_{si}_{gi}_{repeat}.json',row)
        print(net['label'],'fixed',len(fixed),'online',len(online),flush=True)
    write(stage/'fixed_records.json',fixed);write(stage/'online_records.json',online)
    sources=['scripts/320_evaluate_fixed_testset_20261112.py','scripts/arb_polynomial_forcing.py','scripts/225_derivative_split_20261105.py','scripts/151_new_pinn_matched_baselines_20261102.py','scripts/arb_cubic_reconstruction.py','scripts/arb_shaped_remainder.py','scripts/arb_deep_network.py','scripts/arb_reconstruction.py','scripts/arb_certification.py','scripts/arb_problem_families.py','scripts/interval_certification.py']
    summary={'networks':len(nets),'planned_networks':ts['expected'],'missing_candidates':ts['failures'],'fixture_sha256':h,'fixed_records':len(fixed),'online_requests':len(online),'config_sha256':sha(a.config),'source_sha256':{p:sha(ROOT/p) for p in sources},'file_sha256':{p.name:sha(p) for p in sorted(stage.glob('*.json'))},'versions':{'python':platform.python_version(),'python_flint':flint.__version__,'flint':flint.__FLINT_VERSION__},'scope':'Locked protocol, references absent from policy. Timings are same-process requests, not full CLI or independent networks; soft limits can overrun by a completed mesh.'}
    write(stage/'summary.json',summary);os.rename(stage,a.out);print(json.dumps({k:v for k,v in summary.items() if k not in ['file_sha256','source_sha256']}))
if __name__=='__main__':main()
