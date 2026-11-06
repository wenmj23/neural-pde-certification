"""Independent sampled derivative diagnostics with omitted-source controls."""
from fractions import Fraction as F
import hashlib,importlib.util,json,os,tempfile
from itertools import product
from pathlib import Path
from mpmath.libmp import to_rational
from flint import arb
from arb_certification import exact_fraction
from arb_reconstruction import point

ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results/2026-11-06';CONFIG=ROOT/'configs/split_poly_diagnostic_2026-11-06.json';OUT=B/'poly_split_diagnostic'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name,path):
    sp=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m


def main():
    if OUT.exists():raise FileExistsError(OUT)
    c=json.loads(CONFIG.read_text());ic=ROOT/c['fixtures_config'];inputs=json.loads(ic.read_text());prior=ROOT/'results/2026-11-05/independent_poly_split';summary=json.loads((prior/'summary.json').read_text());src=ROOT/'scripts/225_derivative_split_20261105.py';poly=ROOT/'scripts/74_polynomial_forcing_prototype.py';oracle_path=ROOT/'scripts/107_polynomial_derivative_crosscheck_20261101.py'
    if sha(prior/'locked_config.json')!=sha(ic) or sha(src)!=summary['split_source_sha256'] or sha(poly)!=summary['polynomial_source_sha256']:raise ValueError('input/prototype')
    split,oracle=load('split',src),load('oracle',oracle_path);code=poly.read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
    if code.count(needle)!=1:raise ValueError('prototype layout')
    ns={'__file__':str(poly)};exec(compile(code.replace(needle,'OUT=None'),str(poly),'exec'),ns);Poly=ns['PolynomialForcing']
    class Omit(Poly):
        def source_derivative(self,x,j):return arb(0)
    cases={}
    for x in inputs['fixtures']:
        fp=ROOT/x['path']
        if sha(fp)!=x['sha256']:raise ValueError('fixture hash')
        nets=json.loads(fp.read_text())['networks']
        if len(nets)!=12:raise ValueError('network count')
        for n in nets:
            if n['label'] in cases:raise ValueError('duplicate label')
            cases[n['label']]=(n,fp,x['sha256'])
    segments=[tuple(map(F,seg)) for seg in c['segments']];rows=[];faults=[]
    for name,(net,fp,h) in cases.items():
        for (left,right),order in product(segments,c['orders']):
            xs={'left':left,'middle':(left+right)/2,'right':right};reference={}
            for where,x in xs.items():
                a,b=(F(*to_rational(oracle.mp_value(net,x,order,dps)._mpf_)) for dps in c['reference_digits']);reference[where]=(a,b,abs(a-b)<=F(c['relative_stability'])*max(F(1),abs(b)))
            for pieces in c['split_pieces']:
                p=split.SplitProblem(Poly(fp,name,h,net['source_coefficients'],c['precision_bits']),pieces);interval=point(left).lower().union(point(right).upper());ball=p.r_third(interval) if order==3 else p.r_fourth(interval);lo,hi=exact_fraction(ball.lower()),exact_fraction(ball.upper())
                for where,x in xs.items():
                    a,b,stable=reference[where];rows.append({'network':name,'fixture_sha256':h,'segment':[str(left),str(right)],'order':order,'pieces':pieces,'location':where,'x':str(x),'lower':str(lo),'upper':str(hi),'reference80':str(a),'reference110':str(b),'stable':stable,'contained':lo<=b<=hi})
            if net['source_name']==c['fault_source_name']:
                x=xs[c['fault_location']];p=split.SplitProblem(Omit(fp,name,h,net['source_coefficients'],c['precision_bits']),c['fault_pieces']);ball=p.r_third(point(x)) if order==3 else p.r_fourth(point(x));lo,hi=exact_fraction(ball.lower()),exact_fraction(ball.upper());a,b,stable=reference[c['fault_location']];faults.append({'network':name,'fixture_sha256':h,'segment':[str(left),str(right)],'order':order,'x':str(x),'lower':str(lo),'upper':str(hi),'reference80':str(a),'reference110':str(b),'stable':stable,'detected':not lo<=b<=hi})
        print(name,'done',flush=True)
    expected=len(cases)*len(segments)*len(c['orders'])*len(c['split_pieces'])*len(c['points']);fault_n=sum(n['source_name']==c['fault_source_name'] for n,fp,h in cases.values())*len(segments)*len(c['orders'])
    if len(rows)!=expected or len(faults)!=fault_n:raise ValueError('unexpected row count')
    B.mkdir(parents=True,exist_ok=True);temp=Path(tempfile.mkdtemp(prefix='pde_poly_split_diag_'))
    if os.stat(temp).st_dev!=os.stat(B).st_dev:raise OSError('filesystem')
    stage=temp/'bundle';stage.mkdir();(stage/'locked_config.json').write_bytes(CONFIG.read_bytes());(stage/'records.json').write_text(json.dumps(rows,indent=2)+'\n');(stage/'faults.json').write_text(json.dumps(faults,indent=2)+'\n')
    s={'expected':expected,'records':len(rows),'stable':sum(r['stable'] for r in rows),'contained':sum(r['contained'] for r in rows),'fault_expected':fault_n,'fault_records':len(faults),'fault_stable':sum(r['stable'] for r in faults),'fault_detected':sum(r['detected'] for r in faults),'records_sha256':sha(stage/'records.json'),'faults_sha256':sha(stage/'faults.json'),'config_sha256':sha(stage/'locked_config.json'),'fixture_config_sha256':sha(ic),'split_source_sha256':sha(src),'polynomial_source_sha256':sha(poly),'oracle_source_sha256':sha(oracle_path),'script_sha256':sha(Path(__file__)),'scope':c['scope']};(stage/'summary.json').write_text(json.dumps(s,indent=2)+'\n');os.rename(stage,OUT);temp.rmdir();print(json.dumps(s),flush=True)
    return int(s['stable']!=expected or s['contained']!=expected or s['fault_stable']!=fault_n or s['fault_detected']!=fault_n)
if __name__=='__main__':raise SystemExit(main())
