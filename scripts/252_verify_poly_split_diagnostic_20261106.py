"""Read back polynomial derivative diagnostics and recompute point checks."""
from fractions import Fraction as F
import hashlib,json,math
from itertools import product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'results/2026-11-06/poly_split_diagnostic';OUT=ROOT/'results/2026-11-06/poly_split_diagnostic_readback.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(x,s):
    if not x:raise ValueError(s)
def main():
    if OUT.exists():raise FileExistsError(OUT)
    c=json.loads((D/'locked_config.json').read_text());s=json.loads((D/'summary.json').read_text());rows=json.loads((D/'records.json').read_text());faults=json.loads((D/'faults.json').read_text());ic=ROOT/c['fixtures_config'];inputs=json.loads(ic.read_text())
    need(sha(D/'records.json')==s['records_sha256'] and sha(D/'faults.json')==s['faults_sha256'] and sha(D/'locked_config.json')==s['config_sha256']==sha(ROOT/'configs/split_poly_diagnostic_2026-11-06.json'),'output hashes')
    need(sha(ic)==s['fixture_config_sha256'] and sha(ROOT/'scripts/225_derivative_split_20261105.py')==s['split_source_sha256'] and sha(ROOT/'scripts/74_polynomial_forcing_prototype.py')==s['polynomial_source_sha256'] and sha(ROOT/'scripts/107_polynomial_derivative_crosscheck_20261101.py')==s['oracle_source_sha256'] and sha(ROOT/'scripts/251_poly_split_diagnostic_20261106.py')==s['script_sha256'],'source hashes')
    cases={}
    for f in inputs['fixtures']:
        path=ROOT/f['path'];need(sha(path)==f['sha256'],'fixture')
        for n in json.loads(path.read_text())['networks']:
            need(n['label'] not in cases,'duplicate');cases[n['label']]=(n,f['sha256'])
    seg=[tuple(x) for x in c['segments']];keys=[(r['network'],tuple(r['segment']),r['order'],r['pieces'],r['location']) for r in rows];expected=set(product(cases,seg,c['orders'],c['split_pieces'],c['points']));need(len(rows)==len(set(keys))==len(expected)==s['expected']==864 and set(keys)==expected,'keys');idx=dict(zip(keys,rows))
    for r in rows:
        need(r['fixture_sha256']==cases[r['network']][1],'fixture key');a,b=map(F,r['segment']);x={'left':a,'middle':(a+b)/2,'right':b}[r['location']];lo,hi,v80,v110=map(F,(r['lower'],r['upper'],r['reference80'],r['reference110']));stable=abs(v80-v110)<=F(c['relative_stability'])*max(F(1),abs(v110));inside=lo<=v110<=hi
        need(F(r['x'])==x and r['stable']==stable and r['contained']==inside and stable and inside,'point arithmetic')
    faultnames=[n for n,(record,h) in cases.items() if record['source_name']==c['fault_source_name']];fkeys=[(r['network'],tuple(r['segment']),r['order']) for r in faults];expectedfault=set(product(faultnames,seg,c['orders']));need(len(faults)==len(set(fkeys))==len(expectedfault)==s['fault_expected']==48 and set(fkeys)==expectedfault,'fault keys')
    for r in faults:
        n,h=cases[r['network']];a,b=map(F,r['segment']);x=(a+b)/2;order=r['order'];lo,hi,v80,v110=map(F,(r['lower'],r['upper'],r['reference80'],r['reference110']));forcing=sum((F(co)*math.factorial(i)/math.factorial(i-order)*x**(i-order) for i,co in enumerate(n['source_coefficients']) if i>=order),F(0));stable=abs(v80-v110)<=F(c['relative_stability'])*max(F(1),abs(v110));detected=not lo<=v110<=hi;normal=idx[(r['network'],tuple(r['segment']),order,4,'middle')]
        need(r['fixture_sha256']==h and F(r['x'])==x and forcing!=0 and r['stable']==stable and stable and r['detected']==detected and detected and (v80,v110)==(F(normal['reference80']),F(normal['reference110'])),'fault arithmetic')
    comparisons=[]
    for n,segment,order in product(cases,seg,c['orders']):
        r2=idx[(n,segment,order,2,'middle')];r4=idx[(n,segment,order,4,'middle')];w2=F(r2['upper'])-F(r2['lower']);w4=F(r4['upper'])-F(r4['lower']);comparisons.append((w4<w2,w4==w2,w4>w2))
    o={'records':len(rows),'stable':sum(r['stable'] for r in rows),'contained':sum(r['contained'] for r in rows),'fault_records':len(faults),'fault_detected':sum(r['detected'] for r in faults),'derivative_interval_pairs':len(comparisons),'four_piece_narrower':sum(x[0] for x in comparisons),'equal':sum(x[1] for x in comparisons),'wider':sum(x[2] for x in comparisons),'records_sha256':sha(D/'records.json'),'faults_sha256':sha(D/'faults.json'),'scope':'post-exit counts, exact stored arithmetic and numerical-point checks; no continuous interval proof'};OUT.write_text(json.dumps(o,indent=2)+'\n');print(json.dumps(o),flush=True)
if __name__=='__main__':main()
