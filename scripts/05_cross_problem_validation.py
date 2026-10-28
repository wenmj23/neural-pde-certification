"""Cross-problem checks of signed residual certification; ordinary float64 only.

Uncommitted working experiment: archive with the day's summary when requested.
"""
from pathlib import Path
import csv
import json
import platform
import time
import numpy as np
import torch
from torch import nn

OUT = Path(__file__).resolve().parents[1] / 'results/generated/05_cross_problem_validation'


class PINN(nn.Module):
    def __init__(self):
        super().__init__()
        self.hidden = nn.Linear(1, 16)
        self.output = nn.Linear(16, 1)

    def forward(self, x):
        return x*(1-x)*self.output(torch.tanh(self.hidden(x)))


def dx(y, x):
    return torch.autograd.grad(y, x, torch.ones_like(y), create_graph=True)[0]


def residual(model, x, k):
    return -dx(dx(model(x), x), x) - (k*np.pi)**2*torch.sin(k*np.pi*x)


def arrays(model):
    return (model.hidden.weight.detach().numpy().flatten(),
            model.hidden.bias.detach().numpy(),
            model.output.weight.detach().numpy().flatten(),
            model.output.bias.detach().numpy().item())


def fields(x, params, k):
    w,b,a,c=params
    t=np.tanh(np.asarray(x)[...,None]*w+b)
    n=c+(t*a).sum(axis=-1)
    n1=((1-t*t)*(a*w)).sum(axis=-1)
    n2=((-2*t+2*t**3)*(a*w*w)).sum(axis=-1)
    n3=((-2+8*t*t-6*t**4)*(a*w**3)).sum(axis=-1)
    q=x*(1-x); qp=1-2*x
    return (q*n, -(q*n2+2*qp*n1-2*n)-(k*np.pi)**2*np.sin(k*np.pi*x),
            -(q*n3+3*qp*n2-6*n1)-(k*np.pi)**3*np.cos(k*np.pi*x))


def poly(t, order):
    if order==1:return 1-t*t
    if order==2:return -2*t+2*t**3
    return -2+8*t*t-6*t**4


def activation_range(lo, hi, order):
    v1,v2=poly(lo,order),poly(hi,order)
    lower,upper=np.minimum(v1,v2),np.maximum(v1,v2)
    critical={1:[0.],2:[-1/np.sqrt(3),1/np.sqrt(3)],
              3:[-np.sqrt(2/3),0.,np.sqrt(2/3)]}[order]
    for t in critical:
        inside=(lo<=t)&(t<=hi)
        lower=np.where(inside,np.minimum(lower,poly(t,order)),lower)
        upper=np.where(inside,np.maximum(upper,poly(t,order)),upper)
    return lower,upper


def multiply(l1,u1,l2,u2):
    vals=np.array([l1*l2,l1*u2,u1*l2,u1*u2])
    return vals.min(axis=0),vals.max(axis=0)


def cosine_range(left, right, frequency):
    lo,hi=frequency*left,frequency*right
    c1,c2=np.cos(lo),np.cos(hi)
    lower,upper=np.minimum(c1,c2),np.maximum(c1,c2)
    has_max=np.ceil(lo/(2*np.pi))<=np.floor(hi/(2*np.pi))
    has_min=np.ceil((lo-np.pi)/(2*np.pi))<=np.floor((hi-np.pi)/(2*np.pi))
    return np.where(has_min,-1.,lower),np.where(has_max,1.,upper)


def residual_derivative_bounds(nodes, params, k):
    w,b,a,_=params
    left,right=nodes[:-1],nodes[1:]
    z1,z2=left[:,None]*w+b,right[:,None]*w+b
    lo,hi=np.tanh(np.minimum(z1,z2)),np.tanh(np.maximum(z1,z2))
    ranges=[]
    for order in (1,2,3):
        dl,du=activation_range(lo,hi,order)
        coeff=a*w**order
        ranges.append((np.minimum(dl*coeff,du*coeff).sum(axis=1),
                       np.maximum(dl*coeff,du*coeff).sum(axis=1)))
    ql=np.minimum(left*(1-left),right*(1-right))
    qu=np.maximum(left*(1-left),right*(1-right))
    qu=np.where((left<=.5)&(right>=.5),.25,qu)
    t3=multiply(ql,qu,*ranges[2])
    t2=multiply(1-2*right,1-2*left,*ranges[1])
    vl=t3[0]+3*t2[0]-6*ranges[0][1]
    vu=t3[1]+3*t2[1]-6*ranges[0][0]
    # r'=-v'''-f': enclose the forcing derivative too, retaining its sign.
    fl,fu=cosine_range(left,right,k*np.pi)
    fl,fu=fl*(k*np.pi)**3,fu*(k*np.pi)**3
    rl,ru=-vu-fu,-vl-fl
    return np.maximum(abs(rl),abs(ru))


def solve_piecewise_constant(nodes, forcing):
    h=np.diff(nodes); mid=(nodes[:-1]+nodes[1:])/2
    slope=np.sum(forcing*h*(1-mid));value=0.;maximum=0.;trace=[0.]
    for width,f in zip(h,forcing):
        candidates=[0.,width]
        if f!=0 and 0<slope/f<width:candidates.append(slope/f)
        maximum=max(maximum,*[abs(value+slope*t-.5*f*t*t) for t in candidates])
        value+=slope*width-.5*f*width*width;slope-=f*width;trace.append(value)
    assert abs(value)<1e-9, value
    return maximum,np.array(trace)


def evaluate(model,k,seed,step):
    params=arrays(model)
    nodes=np.linspace(0,1,1001);h=np.diff(nodes)
    start=time.perf_counter()
    vn,rn,_=fields(nodes,params,k)
    mi=residual_derivative_bounds(nodes,params,k)
    m=(rn[:-1]+rn[1:])/2;delta=mi*h/2
    Z,zt=solve_piecewise_constant(nodes,m)
    D,dt=solve_piecewise_constant(nodes,delta)
    elapsed=(time.perf_counter()-start)*1000
    lower=max(0.,Z-D);upper=Z+D
    eps_cells=np.maximum(abs(rn[:-1]),abs(rn[1:]))+delta
    unsigned,_=solve_piecewise_constant(nodes,eps_cells)
    # Independent, denser reference evaluation; not used in constructing bounds.
    xe=np.linspace(0,1,10001)
    ve,_,_=fields(xe,params,k)
    reference=np.max(abs(ve-np.sin(k*np.pi*xe)))
    reference_gap=eps_cells.max()*(1/10000)**2/8
    reference_upper=reference+reference_gap
    covers_grid=bool(upper+1e-12>=reference)
    lower_consistent=bool(lower<=reference_upper+1e-12)
    # Implementation cross-checks at independent points and with PyTorch AD.
    test=nodes[:-1,None]+h[:,None]*np.linspace(0,1,11)
    _,rt,rpt=fields(test,params,k)
    assert np.all(abs(rpt)<=mi[:,None]+1e-8)
    assert np.all(abs(rt-m[:,None])<=delta[:,None]+1e-9)
    xa=torch.linspace(0,1,257).reshape(-1,1).requires_grad_(True)
    ra=residual(model,xa,k);rpa=dx(ra,xa)
    _,rnp,rpnp=fields(xa.detach().numpy().flatten(),params,k)
    assert np.allclose(ra.detach().numpy().flatten(),rnp,atol=1e-8,rtol=1e-9)
    assert np.allclose(rpa.detach().numpy().flatten(),rpnp,atol=1e-7,rtol=1e-9)
    assert covers_grid and lower_consistent
    assert np.all(abs(vn[[0,-1]])<1e-14)
    assert np.all(abs(vn-np.sin(k*np.pi*nodes)-zt)<=dt+1e-9)
    return dict(k=k,seed=seed,step=step,error_grid=float(reference),
                reference_gap_bound=float(reference_gap),lower=float(lower),upper=float(upper),
                ratio=float(upper/reference),unsigned_ratio=float(unsigned/reference),
                remainder=float(D),residual_grid_max=float(abs(rn).max()),
                local_M_max=float(mi.max()),certification_ms=float(elapsed),
                covers_reference_grid=covers_grid,lower_consistent=lower_consistent)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    torch.set_default_dtype(torch.float64);torch.set_num_threads(1)
    rows=[];weights={};timings=[]
    for k in (1,5):
        for seed in (7,17,27):
            torch.manual_seed(seed)
            model=PINN();optimizer=torch.optim.Adam(model.parameters(),lr=.01)
            base=torch.linspace(0,1,130)[1:-1].reshape(-1,1)
            rows.append(evaluate(model,k,seed,0))
            start=time.perf_counter()
            for step in range(1,2001):
                x=base.detach().clone().requires_grad_(True)
                optimizer.zero_grad(set_to_none=True)
                loss=(residual(model,x,k)**2).mean()
                loss.backward();optimizer.step()
                if step in (100,500,2000):rows.append(evaluate(model,k,seed,step))
            timings.append(dict(k=k,seed=seed,training_plus_validation_s=time.perf_counter()-start))
            weights[f'k{k}_seed{seed}']={n:p.detach().tolist() for n,p in model.state_dict().items()}
            r=rows[-1]
            print(f'k={k}, seed={seed}: error_grid={r["error_grid"]:.6e}, '
                  f'bound={r["upper"]:.6e}, ratio={r["ratio"]:.4f}, '
                  f'unsigned_ratio={r["unsigned_ratio"]:.2f}',flush=True)
            (OUT/'records.json').write_text(json.dumps(rows,indent=2))
    config=dict(python=platform.python_version(),torch=torch.__version__,dtype='float64',
                domain=[0,1],boundary=[0,0],exact='sin(k*pi*x)',forcing='(k*pi)^2*sin(k*pi*x)',
                k=[1,5],seeds=[7,17,27],width=16,steps=2000,lr=.01,train_points=128,
                residual_check_points=1001,reference_points=10001,checkpoints=[0,100,500,2000],
                caveat='ordinary float64; not outward-rounded; empirical checks do not replace proofs')
    (OUT/'config.json').write_text(json.dumps(config,indent=2))
    (OUT/'final_weights.json').write_text(json.dumps(weights,indent=2))
    (OUT/'timings.json').write_text(json.dumps(timings,indent=2))
    with (OUT/'records.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    summary=dict(cases=len(rows),covered=sum(r['covers_reference_grid'] for r in rows),
                 lower_consistent=sum(r['lower_consistent'] for r in rows),
                 ratio_min=min(r['ratio'] for r in rows),ratio_max=max(r['ratio'] for r in rows))
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
