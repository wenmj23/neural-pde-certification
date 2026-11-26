"""Post-evaluation continuous references on the mpmath interval backend.

Known solutions are used only here, after all fixed-protocol evaluations have
finished. Disjoint or unconverged results are retained and reported.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,importlib.util,json,os,tempfile

ROOT=Path(__file__).resolve().parents[1]
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
Base=module('reference24','scripts/24_depth_continuous_reference.py').MPDeepReference
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(path,obj):
    with path.open('x') as f:json.dump(obj,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def ep(r,side):return F(int(r[side]['numerator']),int(r[side]['denominator']))

class PolyReference(Base):
    def __init__(self,net):
        super().__init__(net)
        self.coefficients=[F(x) for x in net['source_coefficients']]
    def source(self,x):
        c=self.ctx;value=c.mpf(0)
        for coef in reversed(self.coefficients):value=value*x+c.mpf(coef.numerator)/coef.denominator
        return value
    def error(self,x):
        if x not in self.node_cache:
            xx=self.point(x);n,_,_=self.network(xx);c=self.ctx
            u=sum((c.mpf(k.numerator)/k.denominator*(xx-xx**(i+2))/((i+1)*(i+2)) for i,k in enumerate(self.coefficients)),c.mpf(0))
            self.node_cache[x]=xx*(1-xx)*n-u
        return self.node_cache[x]
    def residual(self,left,right):
        self.residual_calls+=1;c=self.ctx;x=c.mpf([self.point(left).a,self.point(right).b]);n,n1,n2=self.network(x,2)
        return -x*(1-x)*n2-2*(1-2*x)*n1+2*n-self.source(x)

def main():
    p=argparse.ArgumentParser();p.add_argument('--training',type=Path,required=True);p.add_argument('--evaluation',type=Path,required=True);p.add_argument('--config',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    c=json.loads(a.config.read_text());train=json.loads((a.training/'summary.json').read_text());es=json.loads((a.evaluation/'summary.json').read_text());fixture=a.training/'frozen_networks.json'
    if sha(fixture)!=train['file_sha256']['frozen_networks.json'] or sha(fixture)!=es['fixture_sha256']:raise ValueError('fixture identity')
    if sha(a.evaluation/'fixed_records.json')!=es['file_sha256']['fixed_records.json'] or sha(a.evaluation/'online_records.json')!=es['file_sha256']['online_records.json']:raise ValueError('evaluation identity')
    nets=json.loads(fixture.read_text())['networks'];fixed=json.loads((a.evaluation/'fixed_records.json').read_text());online=json.loads((a.evaluation/'online_records.json').read_text());refs=[]
    a.out.parent.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix='.fixed-reference-',dir=a.out.parent))
    for i,net in enumerate(nets):
        ref=Base(net) if net['source_kind']=='sine' else PolyReference(net)
        r=ref.maximum(relative=F(c['reference_relative']),absolute=F(c['reference_absolute']),cap=c['reference_cap'])
        row={'network':net['label'],'source_kind':net['source_kind'],**r};refs.append(row);write(stage/f'reference_{i:02d}.json',row)
        print(net['label'],'converged',r['converged'],'splits',r['splits'],flush=True)
    by={r['network']:r for r in refs};checks=[]
    for group,items in [('fixed',fixed),('online',online)]:
        for i,row in enumerate(items):
            q=by[row['network']];lo,hi=ep(q,'lower'),ep(q,'upper')
            b=row if group=='fixed' else row['enclosure']
            if b is None:checks.append({'group':group,'index':i,'network':row['network'],'status':'no_finite_enclosure','reference_converged':q['converged']});continue
            l,u=F(b['lower']),F(b['upper'])
            checks.append({'group':group,'index':i,'network':row['network'],'status':'finite','reference_converged':q['converged'],'overlaps_reference':l<=hi and lo<=u,'contains_reference':l<=lo<=hi<=u})
    write(stage/'references.json',refs);write(stage/'comparison.json',checks)
    summary={'references':len(refs),'reference_converged':sum(r['converged'] for r in refs),'comparisons':len(checks),'finite':sum(z['status']=='finite' for z in checks),'overlap':sum(z.get('overlaps_reference',False) for z in checks),'contains_reference':sum(z.get('contains_reference',False) for z in checks),'fixture_sha256':sha(fixture),'evaluation_fixed_sha256':sha(a.evaluation/'fixed_records.json'),'evaluation_online_sha256':sha(a.evaluation/'online_records.json'),'config_sha256':sha(a.config),'source_sha256':{p:sha(ROOT/p) for p in ['scripts/321_fixed_testset_references_20261112.py','scripts/24_depth_continuous_reference.py','scripts/18_continuous_reference_crosscheck.py','scripts/interval_certification.py']},'file_sha256':{p.name:sha(p) for p in stage.glob('*.json')},'scope':'Reference uses separate mpmath.iv and known exact solutions after evaluations; still conditional on reference backend and implementation. Nonconvergence/disjoint outcomes are reported, not silently discarded.'}
    write(stage/'summary.json',summary);os.rename(stage,a.out);print(json.dumps({k:v for k,v in summary.items() if k not in ('file_sha256','source_sha256')}))
if __name__=='__main__':main()
