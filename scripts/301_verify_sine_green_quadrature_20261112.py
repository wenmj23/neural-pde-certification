"""Post-exit hash, key and exact stored-decimal checks for sine diagnostics."""
from fractions import Fraction
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'results/2026-11-12/sine_green_quadrature_audit'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(c,m):
    if not c:raise AssertionError(m)
def main():
    c=json.loads((D/'locked_config.json').read_text());s=json.loads((D/'summary.json').read_text())
    fixture=ROOT/c['fixture_path'];old=ROOT/c['records_path']
    need(sha(fixture)==c['fixture_sha256']==s['fixture_sha256'],'fixture')
    need(sha(old)==c['records_sha256']==s['old_records_sha256'],'historical records')
    need(sha(D/'locked_config.json')==s['config_sha256'],'config')
    need(sha(ROOT/'scripts/300_sine_green_quadrature_audit_20261112.py')==s['script_sha256'],'producer')
    need(sha(ROOT/'scripts/296_green_quadrature_audit_20261112.py')==s['helper_sha256'],'helper')
    need(sha(D/'records.json')==s['records_sha256'],'records hash')
    nets=json.loads(fixture.read_text())['networks']; labels=[n['label'] for n in nets]
    records=json.loads((D/'records.json').read_text());keys=[(r['network'],r['point']) for r in records]
    need(len(records)==60 and len(set(keys))==60 and set(keys)=={(n,p) for n in labels for p in c['points']},'keys')
    need(len(s['case_sha256'])==60,'case count')
    for ni,n in enumerate(labels):
        for pi,p in enumerate(c['points']):
            name=f'case_{ni:02d}_{pi:02d}.json'
            need(name in s['case_sha256'] and sha(D/name)==s['case_sha256'][name],'case hash')
            need(json.loads((D/name).read_text())==next(r for r in records if (r['network'],r['point'])==(n,p)),'case equality')
    prior={r['network']:r for r in json.loads(old.read_text()) if r['degree']==c['certificate_degree'] and r['intervals']==c['certificate_mesh']}
    diffs=[];cross=[]
    for row in records:
        need([z['digits'] for z in row['results']]==c['digits'],'precision settings')
        for z in row['results']:
            delta=abs(Fraction(z['green'])-Fraction(z['direct_error']));diffs.append(delta)
            need(delta<Fraction(1,10**30) and z['elapsed_seconds']>0,'pointwise consistency')
        crossdiff=abs(Fraction(row['results'][0]['green'])-Fraction(row['results'][1]['green']))
        need(crossdiff<Fraction(1,10**30),'precision agreement');cross.append(crossdiff)
        upper=Fraction(prior[row['network']]['upper'])
        need(Fraction(row['saved_cubic_upper_fraction'])==upper,'prior bound')
        need(row['direct_below_saved_cubic_upper']==(abs(Fraction(row['results'][-1]['direct_error']))<=upper),'upper flag')
    result={'records':len(records),'case_file_checks':len(s['case_sha256']),
            'point_precision_checks':len(diffs),
            'below_saved_upper':sum(row['direct_below_saved_cubic_upper'] for row in records),
            'max_green_direct_difference_as_saved':float(max(diffs)),
            'max_green_cross_precision_difference_as_saved':float(max(cross)),
            'readback_script_sha256':sha(Path(__file__)),
            'scope':'post-exit hashes, manifests, and arithmetic on stored decimal outputs; finite-point numerical check, not an interval proof'}
    p=D/'independent_readback.json'
    if p.exists():raise FileExistsError(p)
    p.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
