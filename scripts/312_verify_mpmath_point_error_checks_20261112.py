"""Post-exit rational and manifest readback for 48 point enclosures."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'results/2026-11-12/mpmath_point_error_checks'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(q,msg):
    if not q:raise AssertionError(msg)
def main():
    c=json.loads((D/'locked_config.json').read_text());s=json.loads((D/'summary.json').read_text())
    inputs={k:ROOT/c[k+'_path'] for k in ['fixture','baseline','spots']}
    for k,p in inputs.items():need(sha(p)==c[k+'_sha256']==s['inputs'][k],k+' hash')
    need(sha(D/'locked_config.json')==s['config_sha256'] and sha(D/'records.json')==s['records_sha256'],'output hashes')
    need(sha(ROOT/'scripts/311_mpmath_point_error_checks_20261112.py')==s['script_sha256'],'producer')
    need(sha(ROOT/'scripts/interval_certification.py')==s['helper_sha256'],'helper')
    old={(r['network'],r['method']):r for r in json.loads(inputs['baseline'].read_text()) if r['intervals']==c['certificate_mesh']}
    spots={(r['network'],r['kind']):r for r in json.loads(inputs['spots'].read_text())}
    rows=json.loads((D/'records.json').read_text());keys=[(r['network'],r['kind']) for r in rows]
    need(len(rows)==48 and len(set(keys))==48 and set(keys)==set(spots),'manifest')
    for r in rows:
        sp=spots[(r['network'],r['kind'])]
        need(F(r['point'])==F(sp['grid_index'],sp['grid_intervals']),'point')
        lo,hi=F(r['signed_error_lower']),F(r['signed_error_upper']);need(lo<=hi,'point interval')
        a=max(F(0),lo if lo>0 else (-hi if hi<0 else F(0)));b=max(abs(lo),abs(hi))
        need(F(r['absolute_error_lower'])==a and F(r['absolute_error_upper'])==b,'absolute interval')
        v=F(sp['error_110']);need(F(r['high_precision_error'])==v,'point reference')
        need(r['contains_high_precision_error']==(lo<=v<=hi),'point containment flag')
        for m in ['signed_quadratic','signed_cubic']:
            U=F(old[(r['network'],m)]['upper']);need(F(r['saved_upper'][m])==U,'prior upper')
            need(r['point_proved_below_saved_upper'][m]==(b<=U),'upper flag')
            need(r['point_proved_above_saved_upper'][m]==(a>U),'contradiction flag')
    out={'records':len(rows),'high_precision_points_contained':sum(r['contains_high_precision_error'] for r in rows),
         'points_below_quadratic_upper':sum(r['point_proved_below_saved_upper']['signed_quadratic'] for r in rows),
         'points_below_cubic_upper':sum(r['point_proved_below_saved_upper']['signed_cubic'] for r in rows),
         'contradictions_quadratic':sum(r['point_proved_above_saved_upper']['signed_quadratic'] for r in rows),
         'contradictions_cubic':sum(r['point_proved_above_saved_upper']['signed_cubic'] for r in rows),
         'readback_script_sha256':sha(Path(__file__)),
         'scope':'conditional interval arithmetic at finitely many selected rational points, with stored arithmetic and hash readback; not continuous-domain verification'}
    p=D/'independent_readback.json'
    if p.exists():raise FileExistsError(p)
    p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
