"""Exact rational boundary-peak controls, not PINNs or new certificates.

For h=3/n^2, r(x)=max(1-x/h,0) is nonnegative. Its exact
zero-boundary solution is piecewise cubic/linear. Check its true maximum
against the global 1/8 upper bound using Fraction only.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,os,tempfile

ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'results/2026-11-12/nonnegative_residual_exact_controls'
NS=(4,8,16,32,64,128,256)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    if DEST.exists():raise FileExistsError(DEST)
    rows=[]
    for n in NS:
        h=F(3,n*n); x=h*(1-F(1,n)); slope=h/2-h*h/6
        def left(t):return slope*t-t*t/2+t*t*t/(6*h)
        def right(t):return h*h*(1-t)/6
        def left_prime(t):return slope-t+t*t/(2*h)
        def right_prime(t):return -h*h/6
        emax=left(x);global_upper=F(1,8)
        checks={'maximum_inside':F(0)<x<h<F(1),
                'zero_boundaries':left(F(0))==right(F(1))==0,
                'continuous_at_h':left(h)==right(h),
                'derivative_continuous_at_h':left_prime(h)==right_prime(h),
                'critical_point':left_prime(x)==0,
                'initial_slope_positive':left_prime(F(0))>0,
                'tail_slope_negative':right_prime(h)<0,
                'upper_bounds_exact_maximum':global_upper>=emax>0}
        if not all(checks.values()):raise ArithmeticError((n,checks))
        rows.append({'n':n,'support_width':str(h),'maximum_location':str(x),
                     'true_max_error':str(emax),'global_residual_upper':str(global_upper),
                     'global_over_true_max':str(global_upper/emax),
                     'global_over_true_max_approx':float(global_upper/emax),'checks':checks})
    DEST.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix='.nonnegative_control_',dir=DEST.parent))
    (stage/'records.json').write_text(json.dumps(rows,indent=2)+'\n')
    summary={'expected':len(NS),'records':len(rows),'script_sha256':sha(Path(__file__)),
             'records_sha256':sha(stage/'records.json'),
             'scope':'constructed nonnegative triangular residuals with exact rational arithmetic; not trained PINNs or execution of the neural certifier'}
    (stage/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    os.rename(stage,DEST);print(json.dumps(summary))
if __name__=='__main__':main()
