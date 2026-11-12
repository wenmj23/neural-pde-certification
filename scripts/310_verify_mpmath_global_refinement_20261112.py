"""Post-exit hash, mesh-path and rational stopping-rule checks."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,math
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'results/2026-11-12/mpmath_global_refinement'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(q,m):
    if not q:raise AssertionError(m)
def main():
    c=json.loads((D/'locked_config.json').read_text());s=json.loads((D/'summary.json').read_text())
    fixture,prior=ROOT/c['fixture_path'],ROOT/c['prior_path']
    need(sha(fixture)==c['fixture_sha256']==s['fixture_sha256'],'fixture')
    need(sha(prior)==c['prior_sha256']==s['prior_sha256'],'prior')
    for name,expected in [('scripts/309_mpmath_global_refinement_20261112.py',s['script_sha256']),
                          ('scripts/307_mpmath_global_residual_audit_20261112.py',s['model_sha256']),
                          ('scripts/interval_certification.py',s['interval_helper_sha256'])]:need(sha(ROOT/name)==expected,name)
    need(sha(D/'locked_config.json')==s['config_sha256'] and sha(D/'records.json')==s['records_sha256'],'outputs')
    rows=json.loads((D/'records.json').read_text());old={(r['network'],r['mesh']):r for r in json.loads(prior.read_text())}
    need(len(rows)==s['records']==len(s['cases']) and {r['network'] for r in rows}==set(c['selected_networks']),'manifest')
    target=F(c['stop_ratio_to_old_global_bound']);summary={}
    for ni,name in enumerate(c['selected_networks']):
        group=[r for r in rows if r['network']==name];meshes=[r['mesh'] for r in group]
        need(meshes==c['meshes'][:len(meshes)] and len(group)>0,'mesh path')
        need(not any(r['goal_met'] for r in group[:-1]),'premature continuation')
        if not group[-1]['goal_met']:need(len(group)==len(c['meshes']),'budget stopping')
        oldU=F(old[(name,1024)]['old_arb_global_upper'])
        for r in group:
            key=f'case_{ni:02d}_{r["mesh"]}.json'
            need(key in s['cases'] and sha(D/key)==s['cases'][key] and json.loads((D/key).read_text())==r,'case')
            upper=F(r['upper']);need(upper==F(r['max_residual_upper'])/8,'upper arithmetic')
            need(F(r['old_arb_global_upper'])==oldU,'old upper')
            need(r['goal_met']==(upper<=target*oldU),'goal flag')
            need(math.isclose(r['upper_over_old_global'],float(upper/oldU),rel_tol=1e-14),'ratio')
            need(0<=r['max_cell']<r['mesh'] and r['elapsed_seconds']>0,'metadata')
        summary[name]={'meshes':meshes,'ratios':[r['upper_over_old_global'] for r in group],
                       'goal_met':group[-1]['goal_met'],'sum_cell_seconds':sum(r['elapsed_seconds'] for r in group)}
    result={'records':len(rows),'networks':len(summary),'groups':summary,
      'readback_script_sha256':sha(Path(__file__)),
      'scope':'post-exit hashes, exact rational arithmetic and stopping-rule readback; mathematical guarantee remains conditional on IV backend and implementation'}
    out=D/'independent_readback.json'
    if out.exists():raise FileExistsError(out)
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
