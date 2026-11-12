"""Post-exit manifest/hash and arithmetic readback; does not recompute the kernel."""
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import math
import statistics

ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'results/2026-11-12/residual_gap_audit'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require(q,msg):
    if not q:raise AssertionError(msg)

def main():
    cfg=json.loads((D/'locked_config.json').read_text())
    s=json.loads((D/'summary.json').read_text())
    fixture=ROOT/cfg['fixture_path'];bp=ROOT/cfg['baseline_path']
    require(sha(fixture)==cfg['fixture_sha256']==s['fixture_sha256'],'fixture hash')
    require(sha(bp)==cfg['baseline_sha256']==s['baseline_sha256'],'baseline hash')
    require(sha(ROOT/'scripts/294_residual_gap_audit_20261112.py')==s['script_sha256'],'producer hash')
    require(sha(D/'locked_config.json')==s['config_sha256'],'configuration hash')
    for name,h in s['result_sha256'].items():require(sha(D/name)==h,'output hash: '+name)
    nets=json.loads(fixture.read_text())['networks']; names={n['label'] for n in nets}
    rows=json.loads((D/'records.json').read_text());spots=json.loads((D/'spot_checks.json').read_text());syn=json.loads((D/'synthetic.json').read_text())
    keys=[(r['network'],r['grid_intervals']) for r in rows]
    require(len(rows)==60 and len(set(keys))==60 and set(keys)=={(n,g) for n in names for g in cfg['grids']},'grid keys')
    skeys=[(r['network'],r['kind']) for r in spots]
    require(len(spots)==48 and len(set(skeys))==48 and set(skeys)=={(n,k) for n in names for k in cfg['high_precision_spots']},'spot keys')
    require(len(syn)==9 and {x['frequency'] for x in syn}==set(cfg['manufactured_frequencies']),'control keys')
    previous={(r['network'],r['method']):r for r in json.loads(bp.read_text()) if r['intervals']==cfg['reference_mesh']}
    for r in rows:
        for method,u in r['saved_64_bounds'].items():
            require(u==float(Fraction(previous[(r['network'],method)]['upper'])),'saved bound mismatch')
        for q in ('grid_max_abs_error','grid_max_abs_residual','green_abs_grid_max'):
            require(math.isfinite(r[q]) and r[q]>0,'finite positive diagnostic')
        require(math.isclose(r['sampled_global_over_sampled_error'],r['grid_max_abs_residual']/8/r['grid_max_abs_error'],rel_tol=1e-14),'ratio arithmetic')
        require(r['sample_error_below_saved_cubic']==(r['grid_max_abs_error']<=r['saved_64_bounds']['signed_cubic']),'sample upper flag')
        require(r['sample_abs_green_below_saved_spatial']==(r['green_abs_grid_max']<=r['saved_64_bounds']['spatial']),'spatial upper flag')
    for z in spots:
        for i,k in enumerate(['error','residual']):
            a,b=Fraction(z[k+'_80']),Fraction(z[k+'_110'])
            require(abs(a-b)<Fraction(1,10**70),'precision stability')
            require(math.isclose(float(abs(a-b)),float(z['precision_abs_differences'][i]),rel_tol=2e-7,abs_tol=0.0),'stored precision difference')
    for z in syn:
        require(math.isclose(z['exact_global_ratio'],(z['frequency']*math.pi)**2/8,rel_tol=1e-14),'sine formula')
    finest=[r for r in rows if r['grid_intervals']==max(cfg['grids'])]
    summary={'records':len(rows),'spot_checks':len(spots),'synthetic_controls':len(syn),
      'both_signs_at_grid_nodes':sum(r['both_signs_on_grid'] for r in finest),
      'sample_error_below_saved_cubic':sum(r['sample_error_below_saved_cubic'] for r in finest),
      'sample_abs_green_below_saved_spatial':sum(r['sample_abs_green_below_saved_spatial'] for r in finest),
      'median_sampled_global_over_error':statistics.median(r['sampled_global_over_sampled_error'] for r in finest),
      'median_sampled_absolute_green_over_error':statistics.median(r['sampled_abs_green_over_sampled_error'] for r in finest),
      'median_sampled_global_over_absolute_green':statistics.median(r['sampled_global_over_sampled_abs_green'] for r in finest),
      'median_green_point_difference':statistics.median(r['green_direct_max_difference'] for r in finest),
      'max_green_point_difference':max(r['green_direct_max_difference'] for r in finest),
      'max_float_residual_discrepancy':max(z['float_residual_difference'] for z in spots),
      'max_float_error_discrepancy':max(z['float_abs_error_difference'] for z in spots),
      'readback_script_sha256':sha(Path(__file__)),
      'scope':'post-exit hash, manifests, stored arithmetic and precision-string comparisons; not an independent continuous-domain proof or rerun of the numerical kernel'}
    path=D/'independent_readback.json'
    if path.exists():raise FileExistsError(path)
    path.write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary))

if __name__=='__main__':main()
