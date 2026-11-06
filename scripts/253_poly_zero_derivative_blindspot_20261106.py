"""Check intentional forcing omission when higher source derivatives vanish.

Constant/linear sources have zero third and fourth derivatives, so a test
limited to these derivatives cannot detect omission from residual itself.
"""
from fractions import Fraction as F
import hashlib,importlib.util,json
from pathlib import Path
from mpmath.libmp import to_rational
from flint import arb
from arb_certification import exact_fraction
from arb_reconstruction import point

ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results/2026-11-06';D=B/'poly_split_diagnostic';OUT=B/'poly_zero_derivative_blindspot.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name,path):
    sp=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
def main():
    if OUT.exists():raise FileExistsError(OUT)
    c=json.loads((D/'locked_config.json').read_text());s=json.loads((D/'summary.json').read_text());ic=ROOT/c['fixtures_config'];inputs=json.loads(ic.read_text());src=ROOT/'scripts/225_derivative_split_20261105.py';poly=ROOT/'scripts/74_polynomial_forcing_prototype.py';oracle_path=ROOT/'scripts/107_polynomial_derivative_crosscheck_20261101.py'
    if sha(D/'records.json')!=s['records_sha256'] or sha(src)!=s['split_source_sha256'] or sha(poly)!=s['polynomial_source_sha256'] or sha(oracle_path)!=s['oracle_source_sha256']:raise ValueError('input/source hash')
    normal={(r['network'],tuple(r['segment']),r['order']):r for r in json.loads((D/'records.json').read_text()) if r['pieces']==4 and r['location']=='middle'}
    split,oracle=load('split',src),load('oracle',oracle_path);code=poly.read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
    if code.count(needle)!=1:raise ValueError('prototype')
    ns={'__file__':str(poly)};exec(compile(code.replace(needle,'OUT=None'),str(poly),'exec'),ns);Poly=ns['PolynomialForcing']
    class Omit(Poly):
        def source_derivative(self,x,j):return arb(0)
    rows=[];nets=0
    for cfg in inputs['fixtures']:
        fp=ROOT/cfg['path']
        if sha(fp)!=cfg['sha256']:raise ValueError('fixture hash')
        for net in json.loads(fp.read_text())['networks']:
            if net['source_name'] not in ('constant_2','linear_6x'):continue
            nets+=1;p=split.SplitProblem(Omit(fp,net['label'],cfg['sha256'],net['source_coefficients'],103),4)
            for raw in c['segments']:
                a,b=map(F,raw);x=(a+b)/2;forcing=sum((F(co)*x**j for j,co in enumerate(net['source_coefficients'])),F(0))
                if not forcing:raise ValueError('source zero at selected point')
                for order in (0,3,4):
                    ball=p.residual(point(x)) if order==0 else p.r_third(point(x)) if order==3 else p.r_fourth(point(x));lo,hi=exact_fraction(ball.lower()),exact_fraction(ball.upper())
                    if order==0:r80,r110=(F(*to_rational(oracle.mp_value(net,x,order,dps)._mpf_)) for dps in (80,110))
                    else:
                        q=normal[(net['label'],tuple(raw),order)];r80,r110=F(q['reference80']),F(q['reference110'])
                    stable=abs(r80-r110)<=F(c['relative_stability'])*max(F(1),abs(r110));detected=not lo<=r110<=hi
                    rows.append({'network':net['label'],'fixture_sha256':cfg['sha256'],'source':net['source_name'],'segment':raw,'x':str(x),'order':order,'source_value':str(forcing),'lower':str(lo),'upper':str(hi),'reference80':str(r80),'reference110':str(r110),'stable':stable,'detected':detected,'expected_detected':order==0})
            print(net['label'],'done',flush=True)
    expected=nets*len(c['segments'])*3
    if nets!=16 or len(rows)!=expected:raise ValueError('request count')
    o={'networks':nets,'expected':expected,'records':len(rows),'stable':sum(x['stable'] for x in rows),'residual_detected':sum(x['detected'] for x in rows if x['order']==0),'high_order_detected':sum(x['detected'] for x in rows if x['order']!=0),'expected_matching':sum(x['detected']==x['expected_detected'] for x in rows),'source_records_sha256':sha(D/'records.json'),'oracle_sha256':sha(oracle_path),'split_sha256':sha(src),'poly_sha256':sha(poly),'script_sha256':sha(Path(__file__)),'rows':rows,'scope':'deliberately omitted source in isolated test; high derivative forcing zero for constant/linear sources; numerical point diagnostic only'}
    tmp=OUT.with_suffix('.tmp')
    with tmp.open('x') as f:json.dump(o,f,indent=2);f.write('\n');f.flush()
    tmp.replace(OUT);print(json.dumps({k:v for k,v in o.items() if k!='rows'}),flush=True)
    return int(o['stable']!=expected or o['expected_matching']!=expected)
if __name__=='__main__':raise SystemExit(main())
