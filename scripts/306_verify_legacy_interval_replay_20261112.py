"""Post-exit readback of an existing mpmath interval-audit replay."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'results/2026-11-12/legacy_mpmath_interval_audit'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ep(r,side):
    a=r[side];x=F(int(a['numerator']),int(a['denominator']));d=F(a['decimal_outward'])
    assert (d<=x if side=='lower' else d>=x)
    return x
def main():
    summary=json.loads((D/'audit_summary.json').read_text());checks=json.loads((D/'checks.json').read_text())
    records=json.loads((D/'certificates.json').read_text());controls=json.loads((D/'analytic_checks.json').read_text())
    ref=json.loads((D/'continuous_reference.json').read_text());dec=json.loads((D/'decisions.json').read_text())
    assert len(checks)==summary['checks_passed']==94 and len({c['name'] for c in checks})==94
    assert all(c['passed'] is True for c in checks)
    assert len(records)==summary['interval_configurations']==6 and len(controls)==summary['analytic_configurations']==13
    assert {(r['precision_digits'],r['intervals']) for r in records}=={(d,n) for d in (30,50) for n in (128,512,1000)}
    for path,h in summary['sources'].items():assert sha(ROOT/path)==h
    case=ROOT/'cases/poisson_1d_sine_k1_seed7'
    assert sha(case/'case.json')==summary['case_sha256'] and sha(case/'weights.json')==summary['weights_sha256']
    rl,ru=ep(ref,'lower'),ep(ref,'upper');assert rl<=ru
    for r in records:assert ep(r,'lower')<=rl<=ru<=ep(r,'upper')
    for r in controls:
        lo,hi=ep(r,'lower'),ep(r,'upper'); name=r['problem']
        if name=='zero': target=F(0)
        elif name=='constant': target=F(1,4)
        elif name=='linear':
            with mp.workdps(80):
                target_mp=2/(3*mp.sqrt(3)); assert mp.mpf(lo.numerator)/lo.denominator<=target_mp<=mp.mpf(hi.numerator)/hi.denominator
            continue
        else:target=F(1,10000)
        assert lo<=target<=hi
    r=next(r for r in records if r['precision_digits']==50 and r['intervals']==512)
    lo,hi=ep(r,'lower'),ep(r,'upper')
    for z in dec:
        tol=F(int(z['tolerance_numerator']),int(z['tolerance_denominator']))
        status='accept' if hi<=tol else ('reject' if lo>tol else 'uncertain')
        assert status==z['status']
    names=['audit_summary.json','checks.json','certificates.json','analytic_checks.json','continuous_reference.json','decisions.json']
    result={'checks':len(checks),'analytic_controls':len(controls),'interval_configurations':len(records),
       'continuous_reference_contained':len(records),'decision_statuses':[z['status'] for z in dec],
       'file_sha256':{n:sha(D/n) for n in names},'readback_script_sha256':sha(Path(__file__)),
       'scope':'post-exit count, input/source hashes, rational endpoint/export and stored-decision checks; one irrational control checked numerically; same interval backend, no independent formal audit'}
    out=D/'independent_readback.json'
    if out.exists():raise FileExistsError(out)
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
