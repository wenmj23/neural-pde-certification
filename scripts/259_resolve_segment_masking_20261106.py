"""Recheck nonzero-source injected faults masked by segment-wide ranges at points."""
from collections import Counter
from fractions import Fraction as F
import hashlib,json
from pathlib import Path
from flint import arb
from arb_certification import exact_fraction
from arb_reconstruction import point
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/2026-11-06';D=BASE/'low_order_interval_diagnostic';OUT=BASE/'low_order_segment_masking_points.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    if OUT.exists():raise FileExistsError(OUT)
    c=json.loads((D/'locked_config.json').read_text());s=json.loads((D/'summary.json').read_text());check=json.loads((BASE/'low_order_interval_diagnostic_readback.json').read_text())
    if sha(D/'records.json')!=s['records_sha256'] or s['records_sha256']!=check['records_sha256']:raise ValueError('records hash')
    poly_path=ROOT/'scripts/74_polynomial_forcing_prototype.py'
    if sha(poly_path)!=s['poly_sha256'] or sha(ROOT/c['fixtures_config'])!=s['fixture_config_sha256']:raise ValueError('source hash')
    code=poly_path.read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
    if code.count(needle)!=1:raise ValueError('prototype')
    ns={'__file__':str(poly_path)};exec(compile(code.replace(needle,'OUT=None'),str(poly_path),'exec'),ns);Poly=ns['PolynomialForcing']
    class Omit(Poly):
        def source_derivative(self,x,j):return arb(0)
    class Reverse(Poly):
        def source_derivative(self,x,j):return -super().source_derivative(x,j)
    nets={}
    for f in json.loads((ROOT/c['fixtures_config']).read_text())['fixtures']:
        path=ROOT/f['path']
        if sha(path)!=f['sha256']:raise ValueError('fixture hash')
        for n in json.loads(path.read_text())['networks']:nets[n['label']]=(n,path,f['sha256'])
    selected=[r for r in json.loads((D/'records.json').read_text()) if r['kind']!='normal' and F(r['source_derivative'])!=0 and r['contained']]
    rows=[]
    for r in selected:
        net,path,digest=nets[r['network']];P=Omit if r['kind']=='omit_source' else Reverse;p=P(path,net['label'],digest,net['source_coefficients'],c['precision_bits']);x=point(F(r['x']));ball=p.residual(x) if r['order']==0 else p.r_second(x);lo,hi=exact_fraction(ball.lower()),exact_fraction(ball.upper());ref=F(r['reference110'])
        rows.append({'network':r['network'],'fixture_sha256':digest,'segment':r['segment'],'location':r['location'],'x':r['x'],'order':r['order'],'kind':r['kind'],'segment_lower':r['lower'],'segment_upper':r['upper'],'point_lower':str(lo),'point_upper':str(hi),'reference110':r['reference110'],'detected':not lo<=ref<=hi})
    counts=Counter((r['kind'],'detected' if r['detected'] else 'contained') for r in rows)
    o={'records':len(rows),'detected':sum(r['detected'] for r in rows),'counts':[{'kind':k[0],'outcome':k[1],'count':v} for k,v in sorted(counts.items())],'source_records_sha256':sha(D/'records.json'),'poly_sha256':sha(poly_path),'script_sha256':sha(Path(__file__)),'rows':rows,'scope':'post-hoc point recheck of selected segment-masked synthetic faults; numerical reference'}
    tmp=OUT.with_suffix('.tmp')
    with tmp.open('x') as f:json.dump(o,f,indent=2);f.write('\n');f.flush()
    tmp.replace(OUT);print(json.dumps({k:v for k,v in o.items() if k!='rows'}))
    return int(len(rows)!=107 or o['detected']!=107)
if __name__=='__main__':raise SystemExit(main())
