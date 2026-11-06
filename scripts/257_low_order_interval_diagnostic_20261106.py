"""Test lower-order segment enclosures and sensitivity to injected source errors."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import factorial
import hashlib, importlib.util, json, os, tempfile
from pathlib import Path
from flint import arb
from mpmath.libmp import to_rational
from arb_certification import exact_fraction
from arb_reconstruction import point

ROOT = Path(__file__).resolve().parents[1]
CFG = ROOT/'configs/low_order_interval_diagnostic_2026-11-06.json'
BASE = ROOT/'results/2026-11-06'
OUT = BASE/'low_order_interval_diagnostic'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    if OUT.exists(): raise FileExistsError(OUT)
    c = json.loads(CFG.read_text()); fixture_config = ROOT/c['fixtures_config']
    prior = json.loads((BASE/'multiorder_fault_diagnostic/summary.json').read_text())
    poly_path = ROOT/'scripts/74_polynomial_forcing_prototype.py'; oracle_path = ROOT/'scripts/107_polynomial_derivative_crosscheck_20261101.py'
    if sha(fixture_config)!=prior['fixture_config_sha256'] or sha(poly_path)!=prior['poly_sha256'] or sha(oracle_path)!=prior['oracle_sha256']: raise ValueError('hash mismatch')
    code = poly_path.read_text().split("fixture=ROOT/'configs/")[0]; needle = "OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
    if code.count(needle)!=1: raise ValueError('prototype')
    ns={'__file__':str(poly_path)};exec(compile(code.replace(needle,'OUT=None'),str(poly_path),'exec'),ns); Poly=ns['PolynomialForcing']
    spec=importlib.util.spec_from_file_location('oracle107',oracle_path);oracle=importlib.util.module_from_spec(spec);spec.loader.exec_module(oracle)
    class Omit(Poly):
        def source_derivative(self,x,j):return arb(0)
    class Reverse(Poly):
        def source_derivative(self,x,j):return -super().source_derivative(x,j)
    rows=[]
    for fixture in json.loads(fixture_config.read_text())['fixtures']:
        path=ROOT/fixture['path']
        if sha(path)!=fixture['sha256']:raise ValueError('fixture digest')
        for net in json.loads(path.read_text())['networks']:
            args=(path,net['label'],fixture['sha256'],net['source_coefficients'],c['precision_bits']);problems={'normal':Poly(*args),'omit_source':Omit(*args),'reverse_source_sign':Reverse(*args)}
            for segment,order in product(c['segments'],c['orders']):
                a,b=map(F,segment); whole=point(a).lower().union(point(b).upper());positions={'left':a,'middle':(a+b)/2,'right':b}; ranges={}
                for kind,p in problems.items():
                    ball=p.residual(whole) if order==0 else p.r_second(whole);ranges[kind]=(exact_fraction(ball.lower()),exact_fraction(ball.upper()))
                for location,x in positions.items():
                    r80,r110=(F(*to_rational(oracle.mp_value(net,x,order,dps)._mpf_)) for dps in c['reference_digits']);stable=abs(r80-r110)<=F(c['relative_stability'])*max(F(1),abs(r110));src=sum((F(v)*factorial(j)/factorial(j-order)*x**(j-order) for j,v in enumerate(net['source_coefficients']) if j>=order),F(0))
                    for kind,(lo,hi) in ranges.items():
                        rows.append({'fixture_sha256':fixture['sha256'],'network':net['label'],'source':net['source_name'],'segment':segment,'location':location,'x':str(x),'order':order,'kind':kind,'source_derivative':str(src),'lower':str(lo),'upper':str(hi),'reference80':str(r80),'reference110':str(r110),'stable':stable,'contained':lo<=r110<=hi})
            print(net['label'],'done',flush=True)
    counts=Counter((r['kind'],'zero' if F(r['source_derivative'])==0 else 'nonzero','contained' if r['contained'] else 'detected') for r in rows)
    stage_parent=Path(tempfile.mkdtemp(prefix='.low_order_',dir=BASE));stage=stage_parent/'bundle';stage.mkdir();(stage/'locked_config.json').write_bytes(CFG.read_bytes());(stage/'records.json').write_text(json.dumps(rows,indent=2)+'\n')
    summary={'records':len(rows),'stable':sum(r['stable'] for r in rows),'normal_contained':sum(r['contained'] for r in rows if r['kind']=='normal'),'counts':[{'kind':k[0],'source_derivative':k[1],'outcome':k[2],'count':v} for k,v in sorted(counts.items())],'records_sha256':sha(stage/'records.json'),'config_sha256':sha(CFG),'fixture_config_sha256':sha(fixture_config),'poly_sha256':sha(poly_path),'oracle_sha256':sha(oracle_path),'script_sha256':sha(Path(__file__)),'scope':c['scope']}
    (stage/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');os.rename(stage,OUT);stage_parent.rmdir();print(json.dumps(summary),flush=True)
    return int(len(rows)!=1296 or summary['stable']!=1296 or summary['normal_contained']!=432)
if __name__=='__main__':raise SystemExit(main())
