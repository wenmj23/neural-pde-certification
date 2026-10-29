"""Follow-up diagnostics and paired timing on two additional frozen seeds."""
import argparse
from datetime import datetime,timezone
from fractions import Fraction
import importlib
import importlib.util
import json
from pathlib import Path
import statistics
import sys

import mpmath as mp
import flint
from flint import arb
from interval_certification import context,endpoint_fraction,sha256
from interval_reformulated import ReformulatedPoisson
from arb_certification import ArbPoisson,exact_fraction

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--study',type=Path);args=parser.parse_args()
    if args.study:study=args.study.resolve()
    else:
        candidates=[p.parent for p in (ROOT/'results/generated/15_interval_optimization').glob('*/config.json')
                    if not (p.parent/'SUPERSEDED.json').exists()]
        study=max(candidates,key=lambda p:p.name)
    out=study/('secondary_followup_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
    out.mkdir(exist_ok=False);checks=[]
    def check(name,condition):
        checks.append(dict(name=name,passed=bool(condition)))
        if not condition:raise AssertionError(name)
    spec=importlib.util.spec_from_file_location('study15',ROOT/'scripts/15_interval_optimization_study.py')
    helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    frac=helper.frac
    main_case=ROOT/'cases/poisson_1d_sine_k1_seed7'
    native=ArbPoisson(main_case,103,'critical');old=ReformulatedPoisson(main_case,context(50),'critical')
    def polynomial(t,j):
        return {2:-2*t+2*t**3,3:-2+8*t**2-6*t**4,4:16*t-40*t**3+24*t**5}[j]
    with mp.workdps(100):
        for a,b in [(Fraction(-1),Fraction(1)),(Fraction(-4,5),Fraction(-1,10)),
                    (Fraction(-1,20),Fraction(1,20)),(Fraction(1,4),Fraction(3,4))]:
            ar=(arb(a.numerator)/a.denominator).union(arb(b.numerator)/b.denominator)
            ml=old.ctx.mpf(a.numerator)/a.denominator;mr=old.ctx.mpf(b.numerator)/b.denominator
            mi=old.ctx.mpf([ml.a,mr.b])
            for j in (2,3,4):
                an=native._range(ar,j);mn=old._range_poly(mi,j)
                for i in range(17):
                    rational=a+(b-a)*i/16;t=mp.mpf(rational.numerator)/rational.denominator
                    v=Fraction(*mp.libmp.to_rational(polynomial(t,j)._mpf_))
                    check(f'poly_range_arb_{a}_{b}_{j}_{i}',exact_fraction(an.lower())<=v<=exact_fraction(an.upper()))
                    check(f'poly_range_mp_{a}_{b}_{j}_{i}',endpoint_fraction(mn)<=v<=endpoint_fraction(mn,1))
        roots={2:[-1/mp.sqrt(3),1/mp.sqrt(3)],3:[-mp.sqrt(mp.mpf(2)/3),mp.mpf(0),mp.sqrt(mp.mpf(2)/3)],
               4:[-mp.sqrt((15+mp.sqrt(105))/30),-mp.sqrt((15-mp.sqrt(105))/30),
                   mp.sqrt((15-mp.sqrt(105))/30),mp.sqrt((15+mp.sqrt(105))/30)]}
        for j,points in roots.items():
            for k,t in enumerate(points):
                n=native._range(native.critical[j][k][0],j)
                m=old._range_poly(old._critical[j][k][0],j)
                v=Fraction(*mp.libmp.to_rational(polynomial(t,j)._mpf_))
                check(f'critical_boundary_arb_{j}_{k}',exact_fraction(n.lower())<=v<=exact_fraction(n.upper()))
                check(f'critical_boundary_mp_{j}_{k}',endpoint_fraction(m)<=v<=endpoint_fraction(m,1))
    print('Polynomial-range and stationary-boundary diagnostics passed.',flush=True)
    reference=json.loads((study/'continuous_references.json').read_text());runs=[]
    methods=[('mp_cached','pilot'),('arb_factored','pilot')];target=Fraction(51,50)
    for repeat in range(5):
        for seed in (17,27):
            case=study/'fixtures'/f'k1_seed{seed}'
            for variant,controller in (methods if repeat%2==0 else list(reversed(methods))):
                r=helper.timed_request(case,variant,controller,target);r.update(seed=seed,repeat=repeat)
                ref=reference[f'k1_seed{seed}']
                check(f'secondary_ref_{seed}_{repeat}_{variant}',frac(r['lower'])<=frac(ref['lower']) and frac(r['upper'])>=frac(ref['upper']))
                check(f'secondary_stop_{seed}_{repeat}_{variant}',frac(r['upper'])<=target*frac(r['lower']))
                runs.append(r);(out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n')
                print('seed',seed,'repeat',repeat+1,variant,'n',r['intervals'],'ms',round(r['wall_ms'],3),flush=True)
    summary=[]
    for seed in (17,27):
        for variant,_ in methods:
            rows=[r for r in runs if r['seed']==seed and r['variant']==variant]
            check(f'secondary_grid_consistent_{seed}_{variant}',len({r['intervals'] for r in rows})==1)
            summary.append(dict(seed=seed,variant=variant,intervals=rows[0]['intervals'],
                median_wall_ms=statistics.median(r['wall_ms'] for r in rows),
                min_wall_ms=min(r['wall_ms'] for r in rows),max_wall_ms=max(r['wall_ms'] for r in rows)))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (out/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    binary=Path(importlib.import_module(arb.__module__).__file__)
    metadata=dict(study=str(study.relative_to(ROOT)),checks_passed=len(checks),timed_requests=len(runs),
        repeats=5,python_flint=flint.__version__,flint=flint.__FLINT_VERSION__,native_extension=binary.name,
        native_extension_sha256=sha256(binary),sources={str(p.relative_to(ROOT)):sha256(p) for p in
        [Path(__file__),ROOT/'scripts/15_interval_optimization_study.py',ROOT/'scripts/arb_certification.py',ROOT/'scripts/interval_reformulated.py']},
        scope='Two more k=1 frozen 16-unit networks, same target 1.02 and controller; no architecture/PDE generalization claim. Sampled high-precision polynomial checks are diagnostics, not a backend proof.')
    (out/'config.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print('FOLLOWUP_DIRECTORY',out,flush=True)


if __name__=='__main__':main()
