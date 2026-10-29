"""Validation-only mpmath interval references for all frozen depth fixtures.

Reuse the independently implemented continuous branch/refine maximum algorithm
from study18; network value/first/second derivatives use direct chain rules.
"""
from datetime import datetime,timezone
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path

from interval_certification import context

ROOT=Path(__file__).resolve().parents[1]
FIXTURE=ROOT/'configs/frozen_networks_depth_2026-10-29.json'
OLD_DRIVER=ROOT/'scripts/18_continuous_reference_crosscheck.py'
spec=importlib.util.spec_from_file_location('reference18',OLD_DRIVER)
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)


class MPDeepReference(base.MPReference):
    def __init__(self,row):
        self.ctx=context(50);c=self.ctx;self.name=row['label'];self.k=row['frequency'];self.omega=self.k*c.pi
        self.layers=[([[c.mpf(float(w)) for w in weights] for weights in layer['weight']],
                      [c.mpf(float(b)) for b in layer['bias']]) for layer in row['layers']]
        self.node_cache={};self.residual_calls=0

    def network(self,x,order=0):
        c=self.ctx;zero=c.mpf(0);one=c.mpf(1);state=[(x,one,zero)]
        for layer,(weights,biases) in enumerate(self.layers):
            nxt=[]
            for w,b in zip(weights,biases):
                z=b+sum((ww*y[0] for ww,y in zip(w,state)),zero)
                z1=sum((ww*y[1] for ww,y in zip(w,state)),zero) if order else zero
                z2=sum((ww*y[2] for ww,y in zip(w,state)),zero) if order else zero
                if layer<len(self.layers)-1:
                    t=1-2/(c.exp(2*z)+1);s=1-t*t
                    nxt.append((t,s*z1,s*z2-2*t*s*z1*z1))
                else:nxt.append((z,z1,z2))
            state=nxt
        return state[0]


def main():
    out=ROOT/'results/generated/24_depth_reference'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True,exist_ok=False)
    refs={};rows=json.loads(FIXTURE.read_text())['networks']
    sources=[Path(__file__),OLD_DRIVER,ROOT/'scripts/interval_certification.py',FIXTURE]
    config=dict(status='running',relative_tolerance='1e-8',absolute_tolerance='1e-14',
                used_for='Validation only; known sin(k*pi*x) never used by certifier',
                scope='Continuous whole-domain enclosure from endpoint errors and residual curvature, different backend from Arb; still conditional on both backends',
                sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    for row in rows:
        r=MPDeepReference(row).maximum(relative=Fraction(1,10**8),absolute=Fraction(1,10**14),cap=5000)
        refs[row['label']]=r
        (out/'references.json').write_text(json.dumps(refs,indent=2)+'\n')
        print('REFERENCE',row['label'],r['lower']['decimal_outward'],r['upper']['decimal_outward'],
              'splits',r['splits'],'seconds',round(r['elapsed_s'],2),'converged',r['converged'],flush=True)
        if not r['converged']:raise RuntimeError('Reference cap exhausted')
    config.update(status='complete',references=len(refs))
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    print('COMPLETED',out,flush=True)


if __name__=='__main__':main()
