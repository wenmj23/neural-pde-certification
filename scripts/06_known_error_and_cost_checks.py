"""Known-error tests and scoped 1D cost comparison; not machine-validated numerics.

Working results are not committed until the requested daily summary.
Requires scipy in addition to the baseline requirements.
"""
from pathlib import Path
import importlib.util
import json
import statistics
import time
import numpy as np
import scipy
from scipy.linalg import solve_banded

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/generated/06_known_error_and_cost_checks'
spec=importlib.util.spec_from_file_location('cross_validation',Path(__file__).with_name('05_cross_problem_validation.py'))
helper=importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)


def known_error_checks():
    # Fixed PDE: -u''=pi^2 sin(pi*x). Candidate v=u*+A sin(k*pi*x).
    # Exact continuous maximum error is A, not a sampled reference value.
    # This tests the signed-residual engine with analytic derivative bounds;
    # network-specific derivative bound generation is checked separately in 05.
    rows=[]
    for amplitude in (1e-2,1e-4,1e-6):
        for k in (1,5,20):
            for intervals in (100,1000,10000):
                x=np.linspace(0,1,intervals+1);h=np.diff(x);omega=k*np.pi
                r=amplitude*omega**2*np.sin(omega*x)
                cl,cu=helper.cosine_range(x[:-1],x[1:],omega)
                mi=amplitude*omega**3*np.maximum(abs(cl),abs(cu))
                m=(r[:-1]+r[1:])/2;delta=mi*h/2
                Z,_=helper.solve_piecewise_constant(x,m)
                D,_=helper.solve_piecewise_constant(x,delta)
                lower=max(0.,Z-D);upper=Z+D
                covered=lower<=amplitude*(1+1e-11) and upper>=amplitude*(1-1e-11)
                assert covered
                rows.append(dict(amplitude=amplitude,k=k,intervals=intervals,
                                 exact_error=amplitude,lower=lower,upper=upper,
                                 ratio=upper/amplitude,covered=bool(covered)))
    return rows


def finite_difference(intervals):
    x=np.linspace(0.,1.,intervals+1)
    n=intervals-1;h=1/intervals
    band=np.zeros((3,n));band[0,1:]=-1;band[1]=2;band[2,:-1]=-1
    rhs=h*h*np.pi**2*np.sin(np.pi*x[1:-1])
    u=np.zeros(intervals+1);u[1:-1]=solve_banded((1,1),band,rhs)
    return x,u


def median_ms(fn,repeats=30):
    for _ in range(3):fn()
    times=[]
    for _ in range(repeats):
        t=time.perf_counter();fn();times.append((time.perf_counter()-t)*1000)
    return statistics.median(times)


def cost_checks():
    weight_path=ROOT/'results/generated/05_cross_problem_validation/final_weights.json'
    weights=json.loads(weight_path.read_text())['k1_seed7']
    params=(np.array(weights['hidden.weight']).flatten(),np.array(weights['hidden.bias']),
            np.array(weights['output.weight']).flatten(),weights['output.bias'][0])
    def certify():
        x=np.linspace(0,1,1001);h=np.diff(x)
        _,r,_=helper.fields(x,params,1)
        mi=helper.residual_derivative_bounds(x,params,1)
        Z,_=helper.solve_piecewise_constant(x,(r[:-1]+r[1:])/2)
        D,_=helper.solve_piecewise_constant(x,mi*h/2)
        return Z+D
    cert_ms=median_ms(certify)
    xref=np.linspace(0,1,10001)
    v,_,_=helper.fields(xref,params,1)
    pinn_error=float(np.max(abs(v-np.sin(np.pi*xref))))
    fd=[]
    for n in (100,1000):
        ms=median_ms(lambda:finite_difference(n))
        x,u=finite_difference(n)
        error=float(np.max(abs(np.interp(xref,x,u)-np.sin(np.pi*xref))))
        fd.append(dict(intervals=n,median_solve_ms=ms,interpolated_reference_error=error))
    return dict(pinn_reference_error=pinn_error,signed_upper=float(certify()),
                certification_median_ms=cert_ms,finite_difference=fd,repeats=30,
                scope='Same machine, CPU, same 1D Poisson PDE. Certification excludes PINN training and reference checks. FD includes grid/system assembly and solve; reference interpolation/error evaluation excluded. Current Python certificate loops versus compiled band solver; not an optimized universal performance comparison.')


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    rows=known_error_checks();cost=cost_checks()
    (OUT/'known_error_results.json').write_text(json.dumps(rows,indent=2))
    (OUT/'cost_comparison.json').write_text(json.dumps(cost,indent=2))
    (OUT/'environment.json').write_text(json.dumps(dict(numpy=np.__version__,scipy=scipy.__version__,
                 precision='float64; no outward rounding'),indent=2))
    print('Known exact-error cases covered:',sum(r['covered'] for r in rows),'/',len(rows))
    for k in (1,5,20):
        selected=[r for r in rows if r['intervals']==1000 and r['k']==k]
        print('k=',k,'1001-point ratio range=',min(r['ratio'] for r in selected),max(r['ratio'] for r in selected))
    print(json.dumps(cost,indent=2))


if __name__=='__main__':main()
