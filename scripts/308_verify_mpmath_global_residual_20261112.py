"""Post-exit input, manifests and rational-arithmetic checks for IV replay."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,math,statistics
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'results/2026-11-12/mpmath_global_residual_audit'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(q,msg):
    if not q:raise AssertionError(msg)
def main():
    c=json.loads((D/'locked_config.json').read_text());s=json.loads((D/'summary.json').read_text())
    paths={k:ROOT/c[k+'_path'] for k in ['fixture','baseline','spots']}
    for k,p in paths.items():need(sha(p)==c[k+'_sha256']==s['inputs'][k],k+' hash')
    need(sha(D/'locked_config.json')==s['config_sha256'],'config')
    need(sha(ROOT/'scripts/307_mpmath_global_residual_audit_20261112.py')==s['script_sha256'],'producer')
    need(sha(ROOT/'scripts/interval_certification.py')==s['interval_helper_sha256'],'interval helper')
    need(sha(D/'records.json')==s['records_sha256'] and sha(D/'point_checks.json')==s['points_sha256'],'results')
    nets=json.loads(paths['fixture'].read_text())['networks'];names=[n['label'] for n in nets]
    old={(r['network'],r['method']):r for r in json.loads(paths['baseline'].read_text()) if r['intervals']==64}
    oldspots={(r['network'],r['kind']):r for r in json.loads(paths['spots'].read_text())}
    rows=json.loads((D/'records.json').read_text());points=json.loads((D/'point_checks.json').read_text())
    keys=[(r['network'],r['mesh']) for r in rows]
    need(len(rows)==36 and len(set(keys))==36 and set(keys)=={(n,m) for n in names for m in c['meshes']},'rows')
    pkeys=[(r['network'],r['kind']) for r in points]
    need(len(points)==48 and len(set(pkeys))==48 and set(pkeys)==set(oldspots),'point keys')
    need(len(s['cases'])==36,'case count')
    for ni,n in enumerate(names):
        for mesh in c['meshes']:
            key=f'case_{ni:02d}_{mesh}.json';need(key in s['cases'] and sha(D/key)==s['cases'][key],'case hash')
            need(json.loads((D/key).read_text())==next(r for r in rows if (r['network'],r['mesh'])==(n,mesh)),'case data')
    for r in rows:
        q=old[(r['network'],'global')];lo,hi=F(q['reference_lower']),F(q['reference_upper'])
        need(F(r['reference_lower'])==lo and F(r['reference_upper'])==hi and F(r['old_arb_global_upper'])==F(q['upper']),'saved bounds')
        need(F(r['max_residual_upper'])/8==F(r['upper']) and F(r['upper'])>=0,'bound arithmetic')
        need(r['contains_saved_reference']==(hi<=F(r['upper'])),'reference flag')
        need(math.isclose(r['upper_over_reference_midpoint'],float(F(r['upper'])/((lo+hi)/2)),rel_tol=1e-14),'ratio')
        need(0<=r['max_cell']<r['mesh'] and r['elapsed_seconds']>0,'cell metadata')
    for p in points:
        source=oldspots[(p['network'],p['kind'])]
        need(F(p['high_precision_residual'])==F(source['residual_110']),'point source')
        need(p['contains_high_precision_point']==(F(p['lower'])<=F(p['high_precision_residual'])<=F(p['upper'])),'point flag')
    summary={'records':len(rows),'point_checks':len(points),'cases':len(s['cases']),
      'saved_reference_contained':sum(r['contains_saved_reference'] for r in rows),
      'high_precision_points_contained':sum(p['contains_high_precision_point'] for p in points),
      'median_upper_over_reference':{str(m):statistics.median(r['upper_over_reference_midpoint'] for r in rows if r['mesh']==m) for m in c['meshes']},
      'readback_script_sha256':sha(Path(__file__)),
      'scope':'post-exit hash/manifest/rational arithmetic; conditional mpmath interval result and finite-point diagnostics, not independent formal audit'}
    out=D/'independent_readback.json'
    if out.exists():raise FileExistsError(out)
    out.write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
if __name__=='__main__':main()
