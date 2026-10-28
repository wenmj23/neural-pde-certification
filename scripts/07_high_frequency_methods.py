"""Compare constant/linear residual reconstruction, adaptivity and vectorization.

This is a working study, not a claim of novelty or machine-rigorous numerics.
"""
from pathlib import Path
import json,csv,time,statistics,importlib.util
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from certification_core import (SineError,BumpError,NeuralResidual,certificate,
    adaptive_nodes,integrate_piecewise_constant,integrate_piecewise_linear)

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/generated/07_high_frequency_methods'


def loop_constant(nodes,forcing):
    h=np.diff(nodes);mid=(nodes[:-1]+nodes[1:])/2
    slope=np.sum(forcing*h*(1-mid));value=0.;maximum=0.
    for width,f in zip(h,forcing):
        ts=[0.,width]
        if f!=0 and 0<slope/f<width:ts.append(slope/f)
        maximum=max(maximum,*[abs(value+slope*t-f*t*t/2) for t in ts])
        value+=slope*width-f*width*width/2;slope-=f*width
    return maximum


def loop_linear(nodes,r):
    h=np.diff(nodes);slope=sum(width*((1-left)*(r0+r1)/2-width*(r0+2*r1)/6)
        for left,width,r0,r1 in zip(nodes[:-1],h,r[:-1],r[1:]))
    value=0.;maximum=0.
    for width,r0,r1 in zip(h,r[:-1],r[1:]):
        a=(r1-r0)/(2*width);ts=[0.,width]
        roots=np.roots([a,r0,-slope]) if a!=0 else ([slope/r0] if r0!=0 else [])
        for root in roots:
            if abs(np.imag(root))<1e-10 and 0<float(np.real(root))<width:ts.append(float(np.real(root)))
        maximum=max(maximum,*[abs(value+slope*t-r0*t*t/2-(r1-r0)*t**3/(6*width)) for t in ts])
        value+=slope*width-width*width*(2*r0+r1)/6;slope-=width*(r0+r1)/2
    return maximum


def checks():
    rng=np.random.default_rng(291)
    for _ in range(20):
        x=np.r_[0.,np.sort(rng.random(22)),1.]
        r=rng.normal(size=len(x));f=r[:-1]
        assert np.isclose(integrate_piecewise_constant(x,f)[0],loop_constant(x,f),atol=1e-12,rtol=1e-10)
        assert np.isclose(integrate_piecewise_linear(x,r)[0],loop_linear(x,r),atol=1e-12,rtol=1e-10)
    x=np.linspace(0,1,71)
    assert np.isclose(integrate_piecewise_linear(x,2*np.ones_like(x))[0],.25,atol=1e-13)
    assert np.isclose(integrate_piecewise_linear(x,x)[0],1/(9*np.sqrt(3)),atol=1e-13)
    for bad in (np.array([.1,1.]),np.array([0.,.5,.5,1.])):
        try:certificate(SineError(1),bad)
        except ValueError:pass
        else:raise AssertionError('Invalid mesh was not rejected')
    print('Vectorized auxiliary-solves checked against independent loop implementations and analytic solutions.',flush=True)


def record(problem,n,method,adaptive=False,reference=None):
    if adaptive:nodes,mesh=adaptive_nodes(problem,n,method)
    else:nodes=np.linspace(0,1,n+1);mesh=dict(mesh_ms=0.,bound_evaluations=0)
    r=certificate(problem,nodes,method)
    exact=reference if reference is not None else problem.exact_error
    if reference is None:
        assert r['lower']<=exact*(1+1e-9)
        assert r['upper']>=exact*(1-1e-9)
    else:
        assert r['upper']>=reference-1e-10
    return dict(problem=problem.name,intervals=n,method=method,adaptive=adaptive,
        reference_error=float(exact),reference_is_exact=reference is None,
        lower=r['lower'],upper=r['upper'],ratio=r['upper']/exact,remainder=r['D'],
        certificate_ms=r['elapsed_ms'],mesh_ms=mesh['mesh_ms'],
        total_ms=r['elapsed_ms']+mesh['mesh_ms'],
        interval_bound_evaluations=n+mesh['bound_evaluations'])


def neural_problems():
    fixture=json.loads((ROOT/'configs/frozen_networks_round07.json').read_text())
    for row in fixture['networks']:
        ws=row['weights']
        params=(np.array(ws['hidden.weight']).flatten(),np.array(ws['hidden.bias']),
                np.array(ws['output.weight']).flatten(),ws['output.bias'][0])
        yield NeuralResidual(params,row['frequency'],row['label']),row['reference_error_grid']


def neural_derivative_checks():
    import torch
    torch.set_default_dtype(torch.float64)
    for problem,_ in neural_problems():
        w,b,a,c=problem.params
        # Separate tensors per interval ensure that local bounds are compared to points in that cell.
        nodes=np.linspace(0,1,1001)
        xx=nodes[:-1,None]+np.diff(nodes)[:,None]*np.linspace(0,1,11)
        x=torch.tensor(xx,requires_grad=True)
        n=c+(torch.tanh(x[...,None]*torch.tensor(w)+torch.tensor(b))*torch.tensor(a)).sum(dim=-1)
        v=x*(1-x)*n
        def deriv(y):return torch.autograd.grad(y,x,torch.ones_like(y),create_graph=True)[0]
        r=-deriv(deriv(v))-(problem.frequency*np.pi)**2*torch.sin(problem.frequency*np.pi*x)
        rp=deriv(r);rpp=deriv(rp)
        for order,rr in [(1,rp),(2,rpp)]:
            bounds=problem.derivative_bound(nodes,order)
            assert np.all(abs(rr.detach().numpy())<=bounds[:,None]*(1+1e-10)+1e-7)
    for problem in [BumpError(.02),BumpError(.005)]:
        nodes=np.linspace(0,1,1001);xx=nodes[:-1,None]+np.diff(nodes)[:,None]*np.linspace(0,1,11)
        x=torch.tensor(xx,requires_grad=True)
        e=4*problem.amplitude*x*(1-x)*torch.exp(-((x-.5)/problem.width)**2)
        def deriv(y):return torch.autograd.grad(y,x,torch.ones_like(y),create_graph=True)[0]
        r=-deriv(deriv(e));rp=deriv(r);rpp=deriv(rp)
        for order,rr in [(1,rp),(2,rpp)]:
            bounds=problem.derivative_bound(nodes,order)
            assert np.all(abs(rr.detach().numpy())<=bounds[:,None]*(1+1e-10)+1e-7)
    print('Local first/second residual derivative enclosures cross-checked with PyTorch AD.',flush=True)


def median_ms(fn):
    for _ in range(3):fn()
    ts=[]
    for _ in range(30):
        t=time.perf_counter();fn();ts.append((time.perf_counter()-t)*1000)
    return statistics.median(ts)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    checks();neural_derivative_checks();rows=[]
    problems=[SineError(k) for k in (1,5,20,50)]+[BumpError(s) for s in (.05,.02,.005)]
    for problem in problems:
        for n in (32,64,128,256,512,1000,2000):
            for method in ('constant','linear'):rows.append(record(problem,n,method))
            rows.append(record(problem,n,'linear',True))
        selected=[r for r in rows if r['problem']==problem.name and r['intervals']==1000]
        print(problem.name,[(r['method'],r['adaptive'],round(r['ratio'],6)) for r in selected],flush=True)
    for amplitude in (1e-2,1e-6,1e-10,-1e-4):
        for k in (1,20):
            p=SineError(k,amplitude)
            for method in ('constant','linear'):
                row=record(p,1000,method);row['problem']+=f'_A{amplitude}';rows.append(row)
    for problem,ref in neural_problems():
        for n in (64,256,1000):
            for method in ('constant','linear'):rows.append(record(problem,n,method,False,ref))
        for n in (64,256,1000):rows.append(record(problem,n,'linear',True,ref))
        print(problem.name,[(r['method'],r['adaptive'],round(r['ratio'],6)) for r in rows
                            if r['problem']==problem.name and r['intervals']==1000],flush=True)
    # Deliberate aliasing: all residual samples vanish mathematically, but the certified remainder must not.
    for method in ('constant','linear'):
        row=record(SineError(200),100,method);row['problem']='alias_sine_200';rows.append(row)
    costs=[]
    for problem,ref in neural_problems():
        if problem.name!='PINN_k1_seed7':continue
        x=np.linspace(0,1,1001)
        for method in ('constant','linear'):
            costs.append(dict(name=method,median_ms=median_ms(lambda:certificate(problem,x,method))))
        values=problem.residual(x);forcing=(values[:-1]+values[1:])/2
        costs.append(dict(name='old_loop_auxiliary_only',median_ms=median_ms(lambda:loop_constant(x,forcing))))
        costs.append(dict(name='vector_auxiliary_only',median_ms=median_ms(lambda:integrate_piecewise_constant(x,forcing))))
    (OUT/'records.json').write_text(json.dumps(rows,indent=2))
    (OUT/'costs.json').write_text(json.dumps(costs,indent=2))
    with (OUT/'records.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    fig,axes=plt.subplots(1,2,figsize=(11,4.3))
    for ax,name in zip(axes,['sine_20','bump_0.005']):
        for method,adaptive,label in [('constant',False,'Constant / uniform'),('linear',False,'Linear / uniform'),('linear',True,'Linear / adaptive')]:
            rs=[r for r in rows if r['problem']==name and r['method']==method and r['adaptive']==adaptive]
            ax.loglog([r['intervals'] for r in rs],[max(1e-12,r['ratio']-1) for r in rs],'-o',label=label)
        ax.set(title=name,xlabel='Final intervals',ylabel='Upper bound / exact error - 1');ax.grid(alpha=.25);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(OUT/'comparison.png',dpi=150);plt.close(fig)
    print('Total method/problem/grid checks:',len(rows),flush=True)
    print('Median cost measurements:',json.dumps(costs),flush=True)


if __name__=='__main__':main()
