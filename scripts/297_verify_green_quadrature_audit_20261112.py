"""Post-exit hash and stored-arithmetic checks for Green point diagnostics."""
from fractions import Fraction
from pathlib import Path
import hashlib
import json

ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'results/2026-11-12/green_quadrature_audit'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(q,msg):
    if not q:raise AssertionError(msg)

def main():
    c=json.loads((D/'locked_config.json').read_text());s=json.loads((D/'summary.json').read_text())
    fixture=ROOT/c['fixture_path'];baseline=ROOT/c['baseline_path']
    need(sha(fixture)==c['fixture_sha256']==s['fixture_sha256'],'fixture')
    need(sha(baseline)==c['baseline_sha256']==s['baseline_sha256'],'baseline')
    need(sha(D/'locked_config.json')==s['config_sha256'],'configuration')
    need(sha(ROOT/'scripts/296_green_quadrature_audit_20261112.py')==s['script_sha256'],'producer')
    need(sha(D/'records.json')==s['records_sha256'],'records hash')
    nets=json.loads(fixture.read_text())['networks'];records=json.loads((D/'records.json').read_text())
    labels=[n['label'] for n in nets]; expected={(n,p) for n in labels for p in c['points']}
    keys=[(r['network'],r['point']) for r in records]
    need(len(records)==len(expected)==84 and len(set(keys))==84 and set(keys)==expected,'record keys')
    cases=s['case_sha256'];need(len(cases)==84,'case hash count')
    for ni,n in enumerate(labels):
        for pi,p in enumerate(c['points']):
            name=f'case_{ni:02d}_{pi:02d}.json'; need(name in cases and sha(D/name)==cases[name],'case hash')
            got=json.loads((D/name).read_text())
            need(got==next(r for r in records if (r['network'],r['point'])==(n,p)),'case equals record')
    old={(r['network'],r['method']):r for r in json.loads(baseline.read_text()) if r['intervals']==c['grid_intervals_for_saved_certificate']}
    deltas=[];cross=[]
    for r in records:
        need([z['digits'] for z in r['results']]==c['digits'],'precision settings')
        d0,d1=r['results']
        for z in [d0,d1]:
            diff=abs(Fraction(z['green'])-Fraction(z['direct_error']))
            need(diff<Fraction(1,10**30),'pointwise consistency')
            need(Fraction(z['difference'])<Fraction(1,10**30),'reported calculation discrepancy')
            need(z['elapsed_seconds']>0,'timing')
            deltas.append(diff)
        crossdiff=abs(Fraction(d0['green'])-Fraction(d1['green']))
        need(crossdiff<Fraction(1,10**30),'cross-precision agreement')
        cross.append(crossdiff)
        u=Fraction(old[(r['network'],'signed_cubic')]['upper'])
        need(Fraction(r['saved_cubic_upper_fraction'])==u,'saved upper')
        need(r['direct_below_saved_cubic_upper']==(abs(Fraction(d1['direct_error']))<=u),'flag')
    summary={'records':len(records),'case_file_checks':len(cases),'point_precision_checks':len(deltas),
       'below_saved_upper':sum(r['direct_below_saved_cubic_upper'] for r in records),
       'max_green_direct_difference_as_saved':float(max(deltas)),
       'max_green_cross_precision_difference_as_saved':float(max(cross)),
       'readback_script_sha256':sha(Path(__file__)),
       'scope':'post-exit input and per-case hash, record keys, stored decimal arithmetic; finite-point numerical diagnostic only'}
    path=D/'independent_readback.json'
    if path.exists():raise FileExistsError(path)
    path.write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary))
if __name__=='__main__':main()
