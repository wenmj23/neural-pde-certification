"""Online baseline then split-2 fallback on previously studied sine models."""
from fractions import Fraction as F
import hashlib,importlib.util,json,os,random,tempfile,time
from itertools import product
from pathlib import Path
from arb_certification import exact_fraction
from arb_cubic_reconstruction import CubicVerifier,SixthOrderDeep
from arb_shaped_remainder import ShapedVerifier
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results/2026-11-06';CONFIG=ROOT/'configs/online_sine_fallback_2026-11-06.json';OUT=B/'online_sine_fallback'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    if OUT.exists():raise FileExistsError(OUT)
    c=json.loads(CONFIG.read_text());prior=ROOT/c['prior'];ic=ROOT/c['input'];old=json.loads((prior/'summary.json').read_text());i=json.loads(ic.read_text());src=ROOT/'scripts/225_derivative_split_20261105.py';f=ROOT/i['fixture'];ref=ROOT/i['reference']
    if sha(prior/'records.json')!=old['records_sha256'] or sha(prior/'locked_config.json')!=sha(ic) or sha(src)!=old['split_source_sha256'] or sha(f)!=i['fixture_sha256'] or sha(ref)!=i['reference_sha256']:raise ValueError('input hash')
    nets={n['label']:n for n in json.loads(f.read_text())['networks']};refs={q['network']:q for q in json.loads(ref.read_text())}
    if len(nets)!=12 or set(nets)!=set(refs) or not all(q['converged'] for q in refs.values()):raise ValueError('manifest')
    sp=importlib.util.spec_from_file_location('split',src);mod=importlib.util.module_from_spec(sp);sp.loader.exec_module(mod);combinations=list(product(nets,c['degrees'],c['caps']));random.Random(c['seed']).shuffle(combinations);rows=[]
    for name,degree,cap in combinations:
        q=refs[name];rl,ru=F(q['lower']),F(q['upper']);row={'network':name,'degree':degree,'cap':cap,'status':'complete','attempts':[]};total=0
        try:
            for pieces in (1,2):
                p=mod.SplitProblem(SixthOrderDeep(f,name,i['fixture_sha256'],i['bits']),pieces);v=ShapedVerifier(p,'moment',subdivision=i['bernstein_depth']) if degree==2 else CubicVerifier(p,'moment',subdivision=i['bernstein_depth'])
                for mesh in c['mesh_sequence']:
                    if mesh>cap:break
                    v.uniform(mesh);start=time.perf_counter();r=v.certify();dt=time.perf_counter()-start;total+=dt;lo,hi=map(exact_fraction,(r['lower'],r['upper']))
                    if not 0<=lo<=rl<=ru<=hi:raise ArithmeticError('reference')
                    goal=lo>0 and hi<=F(c['gap_ratio'])*lo;row['attempts'].append({'pieces':pieces,'mesh':mesh,'seconds':dt,'cumulative_seconds':total,'lower':str(lo),'upper':str(hi),'goal_met':goal,'derivative_calls_cumulative':p.underlying_derivative_calls})
                    if goal:break
                if row['attempts'][-1]['goal_met']:break
            row['fallback_called']=any(x['pieces']==2 for x in row['attempts']);row['goal_met']=row['attempts'][-1]['goal_met']
        except Exception as e:row.update(status='unresolved',error=type(e).__name__+': '+str(e),fallback_called=any(x['pieces']==2 for x in row['attempts']),goal_met=False)
        rows.append(row);print(name,degree,cap,row['status'],row['fallback_called'],row['goal_met'],flush=True)
    if len(rows)!=72 or len({(r['network'],r['degree'],r['cap']) for r in rows})!=72:raise ValueError('request count')
    B.mkdir(parents=True,exist_ok=True);tmp=Path(tempfile.mkdtemp(prefix='pde_online_sine_'))
    if os.stat(tmp).st_dev!=os.stat(B).st_dev:raise OSError('filesystem')
    stage=tmp/'bundle';stage.mkdir();(stage/'locked_config.json').write_bytes(CONFIG.read_bytes());(stage/'records.json').write_text(json.dumps(rows,indent=2)+'\n');s={'expected':72,'records':len(rows),'complete':sum(r['status']=='complete' for r in rows),'unresolved':[r for r in rows if r['status']!='complete'],'attempts':sum(len(r['attempts']) for r in rows),'fallback_calls':sum(r['fallback_called'] for r in rows),'goal_met':sum(r['goal_met'] for r in rows),'records_sha256':sha(stage/'records.json'),'config_sha256':sha(stage/'locked_config.json'),'script_sha256':sha(Path(__file__)),'prior_records_sha256':sha(prior/'records.json'),'split_source_sha256':sha(src),'scope':c['scope'],'timing':'single randomized run, certify-only sum, excludes setup/startup/I-O/reference'};(stage/'summary.json').write_text(json.dumps(s,indent=2)+'\n');os.rename(stage,OUT);tmp.rmdir();print(json.dumps(s),flush=True)
    return int(bool(s['unresolved']))
if __name__=='__main__':raise SystemExit(main())
