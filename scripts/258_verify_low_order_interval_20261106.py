"""Post-exit manifest and stored-rational readback for segment diagnostics."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import factorial
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; DATA=ROOT/'results/2026-11-06/low_order_interval_diagnostic'; OUT=ROOT/'results/2026-11-06/low_order_interval_diagnostic_readback.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    if OUT.exists():raise FileExistsError(OUT)
    c=json.loads((DATA/'locked_config.json').read_text());s=json.loads((DATA/'summary.json').read_text())
    for key,path in {'records_sha256':DATA/'records.json','config_sha256':DATA/'locked_config.json','fixture_config_sha256':ROOT/c['fixtures_config'],'poly_sha256':ROOT/'scripts/74_polynomial_forcing_prototype.py','oracle_sha256':ROOT/'scripts/107_polynomial_derivative_crosscheck_20261101.py','script_sha256':ROOT/'scripts/257_low_order_interval_diagnostic_20261106.py'}.items():
        if sha(path)!=s[key]:raise ValueError('hash '+key)
    networks={}
    for fixture in json.loads((ROOT/c['fixtures_config']).read_text())['fixtures']:
        path=ROOT/fixture['path']
        if sha(path)!=fixture['sha256']:raise ValueError('fixture hash')
        for net in json.loads(path.read_text())['networks']:
            if net['label'] in networks:raise ValueError('duplicate network')
            networks[net['label']]=(net,fixture['sha256'])
    segments=[tuple(seg) for seg in c['segments']];expected=set(product(networks,segments,c['locations'],c['orders'],c['kinds']));seen=set();counts=Counter();ranges={};refs={};groups=Counter()
    for r in json.loads((DATA/'records.json').read_text()):
        key=(r['network'],tuple(r['segment']),r['location'],r['order'],r['kind'])
        if key not in expected or key in seen:raise ValueError('unexpected case')
        seen.add(key);net,digest=networks[key[0]]
        if r['fixture_sha256']!=digest or r['source']!=net['source_name']:raise ValueError('fixture')
        a,b=map(F,key[1]);x={'left':a,'middle':(a+b)/2,'right':b}[key[2]];src=sum((F(v)*factorial(j)/factorial(j-key[3])*x**(j-key[3]) for j,v in enumerate(net['source_coefficients']) if j>=key[3]),F(0))
        if F(r['x'])!=x or F(r['source_derivative'])!=src:raise ValueError('source arithmetic')
        lo,hi,ref80,ref110=map(F,(r['lower'],r['upper'],r['reference80'],r['reference110']));stable=abs(ref80-ref110)<=F(c['relative_stability'])*max(F(1),abs(ref110));inside=lo<=ref110<=hi
        if lo>hi or not stable or not r['stable'] or r['contained']!=inside or (key[4]=='normal' and not inside) or (not src and not inside):raise ValueError('containment')
        rk=(key[0],key[1],key[3],key[4]);pair=(lo,hi)
        if rk in ranges and ranges[rk]!=pair:raise ValueError('segment range changed across locations')
        ranges[rk]=pair; rr=(key[0],key[1],key[2],key[3]);pair=(ref80,ref110)
        if rr in refs and refs[rr]!=pair:raise ValueError('reference mismatch')
        refs[rr]=pair;ck=(key[4],'zero' if not src else 'nonzero','contained' if inside else 'detected');counts[ck]+=1;groups[(r['source'],key[3],*ck)]+=1
    if seen!=expected or len(networks)!=24:raise ValueError('incomplete records')
    declared={(v['kind'],v['source_derivative'],v['outcome']):v['count'] for v in s['counts']}
    normal=sum(v for k,v in counts.items() if k[0]=='normal' and k[2]=='contained')
    if declared!=counts or s['records']!=len(seen) or s['stable']!=len(seen) or s['normal_contained']!=normal or normal!=432:raise ValueError('summary')
    out={'records':len(seen),'normal_contained':normal,'counts':s['counts'],'groups':[{'source':k[0],'order':k[1],'kind':k[2],'source_derivative':k[3],'outcome':k[4],'count':v} for k,v in sorted(groups.items())],'records_sha256':sha(DATA/'records.json'),'scope':'post-exit readback and exact stored-value checks; numerical point diagnostics do not prove interval correctness'}
    with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in out.items() if k!='groups'}))
if __name__=='__main__':main()
