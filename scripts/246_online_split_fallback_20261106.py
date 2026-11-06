"""Online baseline then optional split-2 fallback, fixed mesh caps."""
from fractions import Fraction as F
import hashlib,importlib.util,json,os,random,tempfile,time
from itertools import product
from pathlib import Path
from arb_certification import exact_fraction
from arb_cubic_reconstruction import CubicVerifier
from arb_shaped_remainder import ShapedVerifier
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results/2026-11-06';CONFIG=ROOT/'configs/online_split_fallback_2026-11-06.json';OUT=B/'online_split_fallback'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def val(x):return F(int(x['numerator']),int(x['denominator']))
def main():
    if OUT.exists():raise FileExistsError(OUT)
    c=json.loads(CONFIG.read_text());prior=ROOT/c['prior'];ic=ROOT/c['inputs'];old=json.loads((prior/'summary.json').read_text());inputs=json.loads(ic.read_text());src=ROOT/'scripts/225_derivative_split_20261105.py';poly=ROOT/'scripts/74_polynomial_forcing_prototype.py'
    if sha(prior/'records.json')!=old['records_sha256'] or sha(prior/'locked_config.json')!=sha(ic) or sha(src)!=old['split_source_sha256'] or sha(poly)!=old['polynomial_source_sha256']:raise ValueError('prior/source hash')
    cases={}
    for f in inputs['fixtures']:
        fp,rp=ROOT/f['path'],ROOT/f['reference']
        if sha(fp)!=f['sha256'] or sha(rp)!=f['reference_sha256']:raise ValueError('fixture/reference hash')
        nets=json.loads(fp.read_text())['networks'];refs=json.loads(rp.read_text())
        if len(nets)!=12 or set(n['label'] for n in nets)!=set(refs):raise ValueError('keys')
        for n in nets:
            q=refs[n['label']]
            if n['label'] in cases or not q['converged']:raise ValueError('network/ref')
            cases[n['label']]=(n,fp,f['sha256'],val(q['lower']),val(q['upper']))
    sp=importlib.util.spec_from_file_location('split',src);mod=importlib.util.module_from_spec(sp);sp.loader.exec_module(mod)
    code=poly.read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
    if code.count(needle)!=1:raise ValueError('prototype')
    ns={'__file__':str(poly)};exec(compile(code.replace(needle,'OUT=None'),str(poly),'exec'),ns)
    combinations=list(product(cases,c['degrees'],c['caps']));random.Random(c['seed']).shuffle(combinations);rows=[]
    for name,degree,cap in combinations:
        n,fp,digest,rl,ru=cases[name];row={'network':name,'degree':degree,'cap':cap,'status':'complete','attempts':[]};total=0
        try:
            for pieces in (c['base_pieces'],c['fallback_pieces']):
                p=mod.SplitProblem(ns['PolynomialForcing'](fp,name,digest,n['source_coefficients'],inputs['bits']),pieces);v=ShapedVerifier(p,'moment',subdivision=inputs['bernstein_depth']) if degree==2 else CubicVerifier(p,'moment',subdivision=inputs['bernstein_depth'])
                for mesh in c['mesh_sequence']:
                    if mesh>cap:break
                    v.uniform(mesh);start=time.perf_counter();r=v.certify();dt=time.perf_counter()-start;total+=dt;lo,hi=map(exact_fraction,(r['lower'],r['upper']))
                    if not 0<=lo<=rl<=ru<=hi:raise ArithmeticError('reference not contained')
                    goal=lo>0 and hi<=F(c['gap_ratio'])*lo;row['attempts'].append({'pieces':pieces,'mesh':mesh,'seconds':dt,'cumulative_seconds':total,'lower':str(lo),'upper':str(hi),'goal_met':goal,'derivative_calls_cumulative':p.underlying_derivative_calls})
                    if goal:break
                if row['attempts'][-1]['goal_met']:break
            row['fallback_called']=any(x['pieces']==c['fallback_pieces'] for x in row['attempts']);row['goal_met']=row['attempts'][-1]['goal_met']
        except Exception as e:row.update(status='unresolved',error=type(e).__name__+': '+str(e),fallback_called=any(x['pieces']==c['fallback_pieces'] for x in row['attempts']),goal_met=False)
        rows.append(row);print(name,degree,cap,row['status'],row['fallback_called'],row['goal_met'],flush=True)
    if len(rows)!=144 or len({(r['network'],r['degree'],r['cap']) for r in rows})!=144:raise ValueError('request count')
    B.mkdir(parents=True,exist_ok=True);tmp=Path(tempfile.mkdtemp(prefix='pde_online_fallback_'))
    if os.stat(tmp).st_dev!=os.stat(B).st_dev:raise OSError('different filesystems')
    stage=tmp/'bundle';stage.mkdir();(stage/'locked_config.json').write_bytes(CONFIG.read_bytes());(stage/'records.json').write_text(json.dumps(rows,indent=2)+'\n');s={'expected':144,'records':len(rows),'complete':sum(r['status']=='complete' for r in rows),'unresolved':[r for r in rows if r['status']!='complete'],'attempts':sum(len(r['attempts']) for r in rows),'fallback_calls':sum(r['fallback_called'] for r in rows),'goal_met':sum(r['goal_met'] for r in rows),'records_sha256':sha(stage/'records.json'),'config_sha256':sha(stage/'locked_config.json'),'script_sha256':sha(Path(__file__)),'prior_records_sha256':sha(prior/'records.json'),'split_source_sha256':sha(src),'polynomial_source_sha256':sha(poly),'scope':c['scope'],'timing':c['timing']};(stage/'summary.json').write_text(json.dumps(s,indent=2)+'\n');os.rename(stage,OUT);tmp.rmdir();print(json.dumps(s),flush=True)
    return int(bool(s['unresolved']))
if __name__=='__main__':raise SystemExit(main())
