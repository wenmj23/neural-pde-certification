"""Produce a fixed, unselected depth benchmark; training is only a fixture source."""
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import platform
import time
import numpy as np
import torch
from torch import nn

ROOT=Path(__file__).resolve().parents[1]


class PINN(nn.Module):
    def __init__(self,depth,width=16):
        super().__init__()
        sizes=[1]+[width]*depth+[1]
        self.layers=nn.ModuleList(nn.Linear(a,b) for a,b in zip(sizes[:-1],sizes[1:]))
    def forward(self,x):
        y=x
        for layer in self.layers[:-1]:y=torch.tanh(layer(y))
        return x*(1-x)*self.layers[-1](y)


def dx(y,x):return torch.autograd.grad(y,x,torch.ones_like(y),create_graph=True)[0]
def residual(model,x,k):return -dx(dx(model(x),x),x)-(k*np.pi)**2*torch.sin(k*np.pi*x)


def main():
    torch.set_default_dtype(torch.float64);torch.set_num_threads(1)
    out=ROOT/'results/generated/22_depth_training'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True,exist_ok=False)
    (out/'executed_source.py').write_bytes(Path(__file__).read_bytes())
    rows=[];start=time.perf_counter()
    config=dict(depths=[1,2,3],width=16,frequencies=[1,5],seeds=[7,17,27],steps=2000,
                optimizer='Adam',learning_rate=0.01,collocation_points=128,hard_boundary='x*(1-x)',
                torch=torch.__version__,python=platform.python_version(),dtype='float64',threads=1,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                selection='All 18 final outputs retained; no best-seed or reference-error selection',
                observation='Dense-grid errors are observations only, not continuous guarantees')
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    for depth in (1,2,3):
        for k in (1,5):
            for seed in (7,17,27):
                torch.manual_seed(seed);model=PINN(depth);optimizer=torch.optim.Adam(model.parameters(),lr=.01)
                base=torch.linspace(0,1,130)[1:-1].reshape(-1,1);history=[];t0=time.perf_counter()
                for step in range(1,2001):
                    x=base.detach().clone().requires_grad_(True);optimizer.zero_grad(set_to_none=True)
                    loss=residual(model,x,k).square().mean()
                    if not torch.isfinite(loss):raise ArithmeticError('Nonfinite training loss')
                    loss.backward();optimizer.step()
                    if step in (1,100,500,1000,2000):history.append(dict(step=step,loss_before_update=float(loss.detach())))
                grid=torch.linspace(0,1,20001).reshape(-1,1)
                with torch.no_grad():observed=float((model(grid)-torch.sin(k*np.pi*grid)).abs().max())
                layers=[dict(weight=l.weight.detach().tolist(),bias=l.bias.detach().tolist(),
                             activation='tanh' if i<depth else 'identity') for i,l in enumerate(model.layers)]
                row=dict(label=f'depth{depth}_k{k}_seed{seed}',depth=depth,width=16,frequency=k,seed=seed,
                         parameter_count=sum(p.numel() for p in model.parameters()),layers=layers,
                         training_history=history,training_wall_s=time.perf_counter()-t0,
                         observed_grid_error=observed,observation_points=20001)
                rows.append(row)
                (out/'networks.json').write_text(json.dumps(dict(config=config,networks=rows),indent=2)+'\n')
                print('FROZEN',row['label'],'error_grid',observed,'seconds',round(row['training_wall_s'],2),flush=True)
    config['total_wall_s']=time.perf_counter()-start
    fixture=dict(config=config,networks=rows)
    (out/'networks.json').write_text(json.dumps(fixture,indent=2)+'\n')
    canonical=ROOT/'configs/frozen_networks_depth_2026-10-29.json'
    if canonical.exists():raise FileExistsError('Preserve existing depth fixture')
    canonical.write_bytes((out/'networks.json').read_bytes())
    (out/'completed.json').write_text(json.dumps(dict(networks=len(rows),fixture=str(canonical.relative_to(ROOT)),sha256=hashlib.sha256(canonical.read_bytes()).hexdigest(),elapsed_s=config['total_wall_s']),indent=2)+'\n')
    print('COMPLETED',out,flush=True)


if __name__=='__main__':main()
