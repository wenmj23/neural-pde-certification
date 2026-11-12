"""Explore how direct interval residual bounds change under refinement."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,importlib.util,json,os,tempfile,time
from interval_certification import context,endpoint_fraction
ROOT=Path(__file__).resolve().parents[1]
CFG=ROOT/'configs/mpmath_global_refinement_2026-11-12.json'
DEST=ROOT/'results/2026-11-12/mpmath_global_refinement'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
sp=importlib.util.spec_from_file_location('iv_model',ROOT/'scripts/307_mpmath_global_residual_audit_20261112.py')
helper=importlib.util.module_from_spec(sp);sp.loader.exec_module(helper)

def main():
    if DEST.exists():raise FileExistsError(DEST)
    c=json.loads(CFG.read_text());fixture=ROOT/c['fixture_path'];prior=ROOT/c['prior_path']
    if sha(fixture)!=c['fixture_sha256'] or sha(prior)!=c['prior_sha256']:raise ValueError('input hash')
    nets={n['label']:n for n in json.loads(fixture.read_text())['networks']}
    prev={(r['network'],r['mesh']):r for r in json.loads(prior.read_text())}
    if len(c['selected_networks'])!=4 or not set(c['selected_networks'])<=set(nets):raise ValueError('selected names')
    target=F(c['stop_ratio_to_old_global_bound'])
    DEST.parent.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix='.iv_refinement_',dir=DEST.parent))
    (stage/'locked_config.json').write_bytes(CFG.read_bytes());rows=[];start=time.monotonic()
    for ni,name in enumerate(c['selected_networks']):
        iv=context(c['precision_digits']);rfun,conv=helper.interval_model(nets[name],iv)
        old=F(prev[(name,1024)]['old_arb_global_upper'])
        if F(prev[(name,1024)]['upper'])<=target*old:raise ValueError('prior already achieved target')
        for mesh in c['meshes']:
            max_r=F(0);cell_index=None;t=time.monotonic()
            for i in range(mesh):
                x=iv.mpf([iv.mpf(i)/mesh,iv.mpf(i+1)/mesh]);v=rfun(x)
                lo,hi=endpoint_fraction(v),endpoint_fraction(v,1);b=max(abs(lo),abs(hi))
                if b>max_r:max_r,cell_index=b,i
            upper=max_r/8;achieved=upper<=target*old
            row={'network':name,'mesh':mesh,'max_residual_upper':str(max_r),'upper':str(upper),
                 'old_arb_global_upper':str(old),'upper_over_old_global':float(upper/old),
                 'goal_met':achieved,'max_cell':cell_index,'exact_weight_imports':conv,
                 'elapsed_seconds':time.monotonic()-t}
            (stage/f'case_{ni:02d}_{mesh}.json').write_text(json.dumps(row,indent=2)+'\n')
            rows.append(row);print(name,mesh,'ratio',round(row['upper_over_old_global'],4),
                                   'goal',achieved,'seconds',round(row['elapsed_seconds'],1),flush=True)
            if achieved:break
    (stage/'records.json').write_text(json.dumps(rows,indent=2)+'\n')
    out={'records':len(rows),'networks':len(c['selected_networks']),'config_sha256':sha(CFG),
         'fixture_sha256':sha(fixture),'prior_sha256':sha(prior),'records_sha256':sha(stage/'records.json'),
         'script_sha256':sha(Path(__file__)),'model_sha256':sha(Path(helper.__file__)),
         'interval_helper_sha256':sha(ROOT/'scripts/interval_certification.py'),
         'cases':{p.name:sha(p) for p in sorted(stage.glob('case_*.json'))},
         'elapsed_seconds':time.monotonic()-start,'scope':c['scope']}
    if len(out['cases'])!=len(rows):raise ValueError('case count')
    (stage/'summary.json').write_text(json.dumps(out,indent=2)+'\n');os.rename(stage,DEST)
    print('COMPLETE',json.dumps({k:v for k,v in out.items() if k!='cases'}),flush=True)
if __name__=='__main__':main()
