"""Independent stdlib-only check of selected report facts and elementary identities.

Recomputes from existing records and exact rational arithmetic; it does not
re-run neural-network certification or independently verify interval bounds.
"""
from fractions import Fraction as F
import hashlib,json,statistics
from pathlib import Path
from itertools import product

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/2026-11-05/midterm_evidence_check.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(x,msg):
    if not x:raise ValueError(msg)
def derivative(poly):return [i*poly[i] for i in range(1,len(poly))]
def value(poly,x):return sum((c*x**i for i,c in enumerate(poly)),F(0))


def main():
    if OUT.exists():raise FileExistsError(OUT)
    path=ROOT/'results/2026-11-02/new_pinn_baselines/accuracy/records.json';summ=path.with_name('summary.json');rows=json.loads(path.read_text());s=json.loads(summ.read_text());need(sha(path)==s['records_sha256'],'baseline records hash')
    names={r['network'] for r in rows};methods=('global','spatial','signed_quadratic','signed_cubic');keys=[(r['network'],r['intervals'],r['method']) for r in rows];need(len(names)==12 and len(keys)==len(set(keys))==144 and set(keys)==set(product(names,(16,32,64),methods)),'baseline keys')
    groups={}
    for m in methods:
        g=[r for r in rows if r['method']==m and r['intervals']==64];ratios=[];accepted=0
        for r in g:
            l,u,rl,ru=map(F,(r['lower'],r['upper'],r['reference_lower'],r['reference_upper']))
            need(r['reference_converged'] and l<=rl<=ru<=u and (rl+ru)>0,'baseline reference')
            ratios.append(u/((rl+ru)/2));accepted+=u<=F(1,100000)
            need(r['accepted']['1/100000']==(u<=F(1,100000)),'acceptance flag')
        groups[m]={'networks':len(g),'reference_contained':len(g),'median_upper_over_reference_midpoint':float(statistics.median(ratios)),'accepted_at_1e-5':accepted}
    e=[F(0),F(1,120),F(0),-F(1,12),F(1,8),-F(1,20)]
    # r=x(x-1/2)(x-1); -e''=r, boundary zero, and samples of r vanish.
    r=[F(0),F(1,2),-F(3,2),F(1)]
    need([-c for c in derivative(derivative(e))]==r,'residual polynomial identity')
    need(value(e,F(0))==value(e,F(1))==0 and all(value(r,x)==0 for x in (F(0),F(1,2),F(1))),'boundary / sample nodes')
    need(value(e,F(1,4))==F(5,4096) and derivative(derivative(derivative(r)))==[F(6)],'counterexample')
    # Integral_0^1 G(x,s) ds = x(1-x)/2. Verify symbolically using
    # x(1-x)^2/2 + (1-x)x^2/2 and its maximum 1/8 at x=1/2.
    for x in (F(0),F(1,7),F(1,2),F(4,5),F(1)):
        integral=(1-x)*x*x/2+x*(1-x)*(1-x)/2
        need(integral==x*(1-x)/2 and integral<=F(1,8),'Green integral')
    # Cubic nodal product after y=t-1/2, z=y^2: q=(z-5/36)^2-1/81.
    # For 0<=z<=1/4, q<=1/144 and q>=-1/81, hence max |q|=1/81.
    for z in (F(0),F(1,36),F(5,36),F(1,4)):
        q=(z-F(1,4))*(z-F(1,36))
        need(q==(z-F(5,36))**2-F(1,81) and -F(1,81)<=q<=F(1,144),'cubic polynomial')
    need((F(5,36)-F(1,4))*(F(5,36)-F(1,36))==-F(1,81),'cubic extremum')
    out={'baseline_64':groups,'green_global_factor':'1/8','counterexample_point_error':'5/4096','true_third_residual_derivative':'6','cubic_nodal_product_max':'1/81','cubic_remainder_denominator':1944,'input_records_sha256':sha(path),'script_sha256':sha(Path(__file__)),'scope':'exact algebra and report data readback; not independent numerical-kernel proof'}
    tmp=OUT.with_suffix('.tmp')
    with tmp.open('x') as f:json.dump(out,f,indent=2);f.write('\n');f.flush()
    tmp.replace(OUT);print(json.dumps(out),flush=True)


if __name__=='__main__':main()
