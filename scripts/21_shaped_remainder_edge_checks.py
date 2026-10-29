"""Additional aliasing and cache checks for the new zero-at-nodes envelope."""
from datetime import datetime,timezone
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from arb_certification import exact_fraction,export
from arb_problem_families import SineError,BumpError
from arb_shaped_remainder import ShapedVerifier

ROOT=Path(__file__).resolve().parents[1]


def main():
    out=ROOT/'results/generated/21_shaped_edge_checks'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True,exist_ok=False);checks=[];records=[]
    def check(name,ok):
        checks.append(dict(name=name,passed=bool(ok)))
        (out/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
        if not ok:raise AssertionError(name)
    # All interpolation nodes are zeros of sin(200*pi*x), while ||e||=1e-4.
    for mode in ('constant','full','moment'):
        v=ShapedVerifier(SineError(200),mode);v.uniform(100);r=v.certify()
        sampled_upper=max(exact_fraction(x.abs_upper()) for x in v.node_cache.values())
        check('aliased_samples_near_zero_'+mode,sampled_upper<Fraction(1,10**20))
        check('aliased_error_not_missed_'+mode,exact_fraction(r['upper'])>=Fraction(1,10000))
        check('aliased_lower_zero_'+mode,exact_fraction(r['lower'])==0)
        records.append(dict(case='alias_sine200_N100',mode=mode,lower=export(r['lower'],'down'),upper=export(r['upper'],'up'),sampled_residual_upper=str(sampled_upper)))
    for name,factory in [('sine20',lambda:SineError(20)),('bump',lambda:BumpError('0.005'))]:
        for mode in ('constant','full','moment'):
            variants=[]
            for cache in (False,True):
                v=ShapedVerifier(factory(),mode,cache=cache)
                for n in (16,32,64):v.uniform(n);r=v.certify()
                variants.append((r['lower'],r['upper'],v.residual_evaluations,v.curvature_evaluations))
            a,b=variants
            check('cache_same_bounds_'+name+'_'+mode,a[:2]==b[:2])
            check('cache_saves_calls_'+name+'_'+mode,b[2]<a[2] and b[3]<=a[3])
            records.append(dict(case=name,mode=mode,cache_off_counts=a[2:],cache_on_counts=b[2:]))
    (out/'records.json').write_text(json.dumps(records,indent=2)+'\n')
    sources=[Path(__file__),ROOT/'scripts/arb_shaped_remainder.py',ROOT/'scripts/arb_reconstruction.py',ROOT/'scripts/arb_problem_families.py',ROOT/'scripts/arb_certification.py']
    config=dict(checks_passed=len(checks),sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},scope='Diagnostic edge checks, not formal software verification')
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    print('COMPLETED',out,'checks',len(checks))


if __name__=='__main__':main()
