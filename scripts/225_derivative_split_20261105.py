"""Research-only interval subdivision for residual derivative bounds.

An original entire cell remains unchanged in interpolation.  Its original
interval endpoints are partitioned into adjacent rational endpoints; the
interval hull of all derivative enclosures covers the original interval.
"""
import argparse
from fractions import Fraction as F
import hashlib,json,os,sys,tempfile
from itertools import product
from pathlib import Path
from arb_certification import exact_fraction
from arb_cubic_reconstruction import CubicVerifier,SixthOrderDeep
from arb_reconstruction import point
from arb_shaped_remainder import ShapedVerifier

ROOT=Path(__file__).resolve().parents[1]; CONFIG=ROOT/'configs/derivative_split_2026-11-05.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(lo,hi,t):return 'within' if hi<=t else 'exceeds' if lo>t else 'unknown'


class SplitProblem:
    def __init__(self,problem,pieces):self.problem=problem;self.pieces=pieces;self.underlying_derivative_calls=0
    def __getattr__(self,key):return getattr(self.problem,key)
    def split_eval(self,x,method):
        f=getattr(self.problem,method)
        if self.pieces==1:
            self.underlying_derivative_calls+=1;return f(x)
        a,b=exact_fraction(x.lower()),exact_fraction(x.upper())
        if a>b:raise ArithmeticError('reversed cell')
        points=[a+(b-a)*F(i,self.pieces) for i in range(self.pieces+1)]
        out=None
        for left,right in zip(points[:-1],points[1:]):
            cell=point(left).lower().union(point(right).upper())
            if exact_fraction(cell.lower())>left or exact_fraction(cell.upper())<right:raise ArithmeticError('partition coverage')
            y=f(cell);self.underlying_derivative_calls+=1
            if not y.is_finite():raise ArithmeticError('nonfinite derivative')
            out=y if out is None else out.union(y)
        return out
    def r_third(self,x):return self.split_eval(x,'r_third')
    def r_fourth(self,x):return self.split_eval(x,'r_fourth')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);out=ap.parse_args().out.resolve()
    out.parent.mkdir(parents=True,exist_ok=True)
    if out.exists():raise FileExistsError(out)
    c=json.loads(CONFIG.read_text());fixture=ROOT/c['fixture'];ref=ROOT/c['reference']
    if sha(fixture)!=c['fixture_sha256'] or sha(ref)!=c['reference_sha256']:raise ValueError('fixed input hash mismatch')
    nets=json.loads(fixture.read_text())['networks'];expected=set(product(c['networks']['frequency'],c['networks']['depth'],c['networks']['seed']))
    if len(nets)!=12 or {(n['frequency'],n['depth'],n['seed']) for n in nets}!=expected:raise ValueError('network manifest')
    temp=Path(tempfile.mkdtemp(prefix='pde_derivative_split_'))
    if os.stat(temp).st_dev!=os.stat(out.parent).st_dev:raise OSError('different filesystem')
    stage=temp/'bundle';stage.mkdir();(stage/'locked_config.json').write_bytes(CONFIG.read_bytes());rows=[]
    for net in nets:
        for deg,mesh,pieces in product(c['degree'],c['mesh'],c['derivative_splits']):
            row={'network':net['label'],'frequency':net['frequency'],'degree':deg,'mesh':mesh,'pieces':pieces}
            try:
                prob=SplitProblem(SixthOrderDeep(fixture,net['label'],c['fixture_sha256'],c['precision_bits']),pieces)
                v=ShapedVerifier(prob,'moment',subdivision=c['bernstein_depth']) if deg==2 else CubicVerifier(prob,'moment',subdivision=c['bernstein_depth']);v.uniform(mesh)
                r=v.certify();lo,hi,zl,zu,du=map(exact_fraction,(r['lower'],r['upper'],r['z_lower'],r['z_upper'],r['remainder_upper']))
                if not (0<=lo<=hi and 0<=zl<=zu and du>=0 and lo<=max(F(0),zl-du) and hi>=zu+du):raise ArithmeticError('invalid enclosure')
                row.update(status='finite',lower=str(lo),upper=str(hi),z_lower=str(zl),z_upper=str(zu),remainder_upper=str(du),gap_met=lo>0 and hi<=F(c['gap_ratio'])*lo,decisions={t:dec(lo,hi,F(t)) for t in c['tolerances']},interval_calls=r['curvature_evaluations'],underlying_derivative_calls=prob.underlying_derivative_calls)
            except Exception as e:row.update(status='unresolved',error=type(e).__name__+': '+str(e))
            rows.append(row)
        print(net['label'],'done',flush=True)
    refs={r['network']:r for r in json.loads(ref.read_text())}
    if set(refs)!={n['label'] for n in nets} or sha(ref)!=c['reference_sha256']:raise ValueError('reference manifest')
    for row in rows:
        if row['status']=='finite':
            q=refs[row['network']];row['reference_converged']=q['converged'];row['reference_contained']=F(row['lower'])<=F(q['lower'])<=F(q['upper'])<=F(row['upper'])
    expected_count=len(nets)*len(c['degree'])*len(c['mesh'])*len(c['derivative_splits'])
    if len(rows)!=expected_count or len({(r['network'],r['degree'],r['mesh'],r['pieces']) for r in rows})!=expected_count:raise ValueError('missing or duplicate')
    (stage/'records.json').write_text(json.dumps(rows,indent=2)+'\n');s={'records':len(rows),'expected':expected_count,'finite':sum(r['status']=='finite' for r in rows),'reference_contained':sum(r.get('reference_contained',False) for r in rows),'unresolved':[r for r in rows if r['status']!='finite'],'records_sha256':sha(stage/'records.json'),'config_sha256':sha(stage/'locked_config.json'),'script_sha256':sha(Path(__file__)),'scope':c['scope']}
    (stage/'summary.json').write_text(json.dumps(s,indent=2)+'\n');os.rename(stage,out);temp.rmdir();print(json.dumps(s),flush=True)
    return int(s['finite']!=expected_count or s['reference_contained']!=expected_count)


if __name__=='__main__':sys.exit(main())
