"""Post-hoc subdivision of previously masked, synthetic source-fault cases."""
from collections import Counter
from fractions import Fraction as F
import hashlib,json
from pathlib import Path
from flint import arb
from arb_certification import exact_fraction
from arb_reconstruction import point
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/2026-11-06';SRC=BASE/'low_order_segment_masking_points.json';OUT=ROOT/'results/2026-11-07/subdivision_fault_masking.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    if OUT.exists():raise FileExistsError(OUT)
    data=json.loads(SRC.read_text());prior=json.loads((BASE/'low_order_segment_masking_points_readback.json').read_text())
    if sha(SRC)!=prior['records_sha256']:raise ValueError('source digest')
    cfg=json.loads((BASE/'low_order_interval_diagnostic/locked_config.json').read_text());pc=ROOT/'scripts/74_polynomial_forcing_prototype.py'
    if sha(pc)!=data['poly_sha256']:raise ValueError('code digest')
    code=pc.read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
    if code.count(needle)!=1:raise ValueError('prototype')
    ns={'__file__':str(pc)};exec(compile(code.replace(needle,'OUT=None'),str(pc),'exec'),ns);Poly=ns['PolynomialForcing']
    class Omit(Poly):
        def source_derivative(self,x,j):return arb(0)
    class Reverse(Poly):
        def source_derivative(self,x,j):return -super().source_derivative(x,j)
    nets={}
    for f in json.loads((ROOT/cfg['fixtures_config']).read_text())['fixtures']:
        p=ROOT/f['path']
        if sha(p)!=f['sha256']:raise ValueError('fixture digest')
        for n in json.loads(p.read_text())['networks']:nets[n['label']]=(n,p,f['sha256'])
    rows=[]
    for r in data['rows']:
        n,p,h=nets[r['network']];args=(p,n['label'],h,n['source_coefficients'],103);bad=(Omit if r['kind']=='omit_source' else Reverse)(*args);good=Poly(*args);a,b=map(F,r['segment']);x=F(r['x']);ref=F(r['reference110'])
        for pieces in (2,4,8,16):
            width=(b-a)/pieces;idx=min((x-a)//width,pieces-1);left=a+idx*width;right=left+width
            if not left<=x<=right:raise ValueError('subinterval')
            domain=point(left).lower().union(point(right).upper())
            vals=[]
            for prob in (good,bad):
                ball=prob.residual(domain) if r['order']==0 else prob.r_second(domain);vals.append((exact_fraction(ball.lower()),exact_fraction(ball.upper())))
            (gl,gh),(bl,bh)=vals
            rows.append({'network':r['network'],'segment':r['segment'],'location':r['location'],'x':str(x),'order':r['order'],'kind':r['kind'],'pieces':pieces,'left':str(left),'right':str(right),'reference':str(ref),'good_lower':str(gl),'good_upper':str(gh),'bad_lower':str(bl),'bad_upper':str(bh),'good_contains':gl<=ref<=gh,'fault_detected':not bl<=ref<=bh})
    counts=Counter((r['pieces'],r['kind'],r['fault_detected']) for r in rows)
    o={'cases':len(data['rows']),'records':len(rows),'good_contained':sum(r['good_contains'] for r in rows),'counts':[{'pieces':k[0],'kind':k[1],'detected':k[2],'count':v} for k,v in sorted(counts.items())],'input_sha256':sha(SRC),'prototype_sha256':sha(pc),'script_sha256':sha(Path(__file__)),'scope':'post-hoc subdivision of 107 selected previously masked synthetic faults; point diagnostics only','rows':rows}
    OUT.parent.mkdir(parents=True,exist_ok=True);tmp=OUT.with_suffix('.tmp')
    with tmp.open('x') as f:json.dump(o,f,indent=2);f.write('\n')
    tmp.replace(OUT);print(json.dumps({k:v for k,v in o.items() if k!='rows'}));return int(len(rows)!=428 or o['good_contained']!=428)
if __name__=='__main__':raise SystemExit(main())
