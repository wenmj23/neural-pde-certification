"""Validation-only continuous maximum enclosure using an independent backend.

Known sin(k*pi*x) is used here, never by certificate construction or stopping.
Branch/refine endpoint error + residual-curvature bounds until the global max
bracket is narrow. Covers the whole domain, including discarded cells. Uses
mpmath interval arithmetic, not the Arb higher-order reconstruction backend.
"""
import argparse
from datetime import datetime,timezone
from fractions import Fraction
import hashlib
import heapq
import json
from pathlib import Path
import sys
import time

from interval_certification import context,endpoint_fraction,export_bounds

ROOT=Path(__file__).resolve().parents[1]
FIXTURE=ROOT/'configs/frozen_networks_round07.json'


def frac(e):return Fraction(int(e['numerator']),int(e['denominator']))


class MPReference:
    def __init__(self,row):
        self.ctx=context(50);c=self.ctx;w=row['weights'];self.name=row['label'];self.k=row['frequency']
        self.w=[c.mpf(float(v[0])) for v in w['hidden.weight']];self.b=[c.mpf(float(v)) for v in w['hidden.bias']]
        self.a=[c.mpf(float(v)) for v in w['output.weight'][0]];self.c=c.mpf(float(w['output.bias'][0]));self.omega=self.k*c.pi
        self.node_cache={};self.residual_calls=0
    def point(self,x):return self.ctx.mpf(x.numerator)/x.denominator
    def network(self,x,order=0):
        n=self.c;n1=self.ctx.mpf(0);n2=self.ctx.mpf(0)
        for w,b,a in zip(self.w,self.b,self.a):
            t=1-2/(self.ctx.exp(2*(w*x+b))+1);n+=a*t
            if order:
                n1+=a*w*(1-t*t);n2+=a*w*w*(-2*t+2*t*t*t)
        return n,n1,n2
    def error(self,x):
        if x not in self.node_cache:
            xx=self.point(x);n,_,_=self.network(xx)
            self.node_cache[x]=xx*(1-xx)*n-self.ctx.sin(self.omega*xx)
        return self.node_cache[x]
    def residual(self,left,right):
        self.residual_calls+=1;c=self.ctx
        x=c.mpf([self.point(left).a,self.point(right).b]);n,n1,n2=self.network(x,2)
        return -x*(1-x)*n2-2*(1-2*x)*n1+2*n-self.omega**2*c.sin(self.omega*x)
    def maximum(self,relative=Fraction(1,10**10),absolute=Fraction(1,10**16),cap=20000):
        start=time.perf_counter();c=self.ctx;lower=c.mpf(0);heap=[];counter=0;splits=0
        def cell(a,b):
            nonlocal lower,counter
            ea,eb=self.error(a),self.error(b)
            lower=max(lower,abs(ea).a,abs(eb).a)
            width=self.point(b)-self.point(a)
            curvature=abs(self.residual(a,b)).b
            upper=(max(abs(ea).b,abs(eb).b)+curvature*width**2/8).b
            heapq.heappush(heap,(-endpoint_fraction(upper),counter,a,b,upper));counter+=1
        for i in range(64):cell(Fraction(i,64),Fraction(i+1,64))
        converged=False
        while heap:
            top=heap[0];uf=-top[0];lf=endpoint_fraction(lower)
            tolerance=max(absolute,relative*lf)
            if uf-lf<=tolerance:converged=True;break
            if splits>=cap:break
            _,_,a,b,_=heapq.heappop(heap)
            # A cell already below the global lower cannot contain the maximum.
            if uf<=lf:continue
            midpoint=(a+b)/2;cell(a,midpoint);cell(midpoint,b);splits+=1
        upper=max(lower,heap[0][4]) if heap else lower
        result=dict(lower=lower,upper=upper)
        return dict(**export_bounds(result),converged=converged,splits=splits,
            error_point_evaluations=len(self.node_cache),residual_interval_evaluations=self.residual_calls,
            elapsed_s=time.perf_counter()-start,bits=self.ctx.prec,domain_coverage='64 initial cells, all active upper bounds retained; discarded cells below global lower',
            known_solution_used=True,used_by_certifier=False)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--study',type=Path);args=parser.parse_args()
    if args.study:study=args.study.resolve()
    else:study=max([p.parent for p in (ROOT/'results/generated/17_reconstruction_adaptivity').glob('*/config.json')],key=lambda p:p.name)
    out=study/('independent_reference_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
    out.mkdir(exist_ok=False);fixture=json.loads(FIXTURE.read_text());refs={};checks=[]
    for row in fixture['networks']:
        result=MPReference(row).maximum()
        if not result['converged']:raise RuntimeError('Continuous reference did not converge: '+row['label'])
        refs[row['label']]=result
        (out/'references.json').write_text(json.dumps(refs,indent=2)+'\n')
        print(row['label'],result['lower']['decimal_outward'],result['upper']['decimal_outward'],
              'splits',result['splits'],'seconds',round(result['elapsed_s'],3),flush=True)
    for filename in ['equal_sample_screen.json','runs.json','extra_runs.json']:
        rows=json.loads((study/filename).read_text())
        for i,r in enumerate(rows):
            if r['problem'] not in refs:continue
            ref=refs[r['problem']]
            overlap=frac(r['lower'])<=frac(ref['upper']) and frac(r['upper'])>=frac(ref['lower'])
            contained=frac(r['lower'])<=frac(ref['lower']) and frac(r['upper'])>=frac(ref['upper'])
            checks.append(dict(file=filename,index=i,problem=r['problem'],overlap=overlap,contains_reference=contained))
            if not overlap:raise AssertionError('Disjoint independent reference and certificate')
    (out/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    meta=dict(scope='Validation-only continuous maximum enclosure, different backend and different formula; both backends still trusted',
        references=len(refs),checks=len(checks),all_overlap=all(r['overlap'] for r in checks),
        all_contain_reference=all(r['contains_reference'] for r in checks),
        contained_count=sum(r['contains_reference'] for r in checks),fixture_sha256=hashlib.sha256(FIXTURE.read_bytes()).hexdigest(),
        sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'scripts/interval_certification.py']})
    (out/'config.json').write_text(json.dumps(meta,indent=2)+'\n')
    print('CROSSCHECK',json.dumps(meta),flush=True)
    print('REFERENCE_DIRECTORY',out,flush=True)


if __name__=='__main__':main()
