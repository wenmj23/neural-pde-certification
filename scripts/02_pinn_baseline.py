# -*- coding: utf-8 -*-
"""Experiment 02: extracted from the executed notebook. See README for assumptions."""
from pathlib import Path
import json
import matplotlib
matplotlib.use("Agg")
OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "generated" / "02_pinn_baseline"
OUT_DIR.mkdir(parents=True, exist_ok=True)

import platform
import time
import numpy as np
import torch
from torch import nn
import matplotlib.pyplot as plt

torch.set_default_dtype(torch.float64)
torch.set_num_threads(1)
torch.manual_seed(7)
np.random.seed(7)
print("Python:", platform.python_version())
print("PyTorch:", torch.__version__)
print("Device: CPU; dtype: float64; seed: 7")


class PoissonPINN(nn.Module):
    def __init__(self, width=16):
        super().__init__()
        self.hidden = nn.Linear(1, width)
        self.output = nn.Linear(width, 1)

    def forward(self, x):
        n = self.output(torch.tanh(self.hidden(x)))
        return x * (1.0 - x) * n

model = PoissonPINN()
endpoints = torch.tensor([[0.0], [1.0]])
print("Number of parameters:", sum(p.numel() for p in model.parameters()))
print("Boundary outputs:", model(endpoints).detach().flatten().tolist())
assert torch.equal(model(endpoints).detach(), torch.zeros_like(endpoints))


def dx(y, x):
    return torch.autograd.grad(
        y, x, grad_outputs=torch.ones_like(y), create_graph=True
    )[0]

def residual(model, x):
    v = model(x)
    v_x = dx(v, x)
    v_xx = dx(v_x, x)
    return -v_xx - 2.0

# 用解析解检查求导与残差的符号；这不是训练数据。
check_x = torch.linspace(0, 1, 17).reshape(-1, 1).requires_grad_(True)
check_u = check_x * (1-check_x)
check_r = -dx(dx(check_u, check_x), check_x) - 2.0
print("Residual of the analytic solution:", check_r.detach().abs().max().item())
assert check_r.detach().abs().max().item() == 0.0


def derivative_bound(model):
    w = model.hidden.weight.detach().reshape(-1)
    a = model.output.weight.detach().reshape(-1)
    b1 = torch.sum(torch.abs(a*w)).item()
    b2 = (2*torch.sum(torch.abs(a)*w.abs()**2)).item()
    b3 = (2*torch.sum(torch.abs(a)*w.abs()**3)).item()
    return 0.25*b3 + 3*b2 + 6*b1

def evaluate(model, step):
    xr = torch.linspace(0, 1, 1001).reshape(-1, 1).requires_grad_(True)
    rr = residual(model, xr)
    r_grid = rr.detach().abs().max().item()
    xe = torch.linspace(0, 1, 10001).reshape(-1, 1)
    with torch.no_grad():
        e_grid = (model(xe) - xe*(1-xe)).abs().max().item()
        bc = model(endpoints).abs().max().item()
    m = derivative_bound(model)
    h = 1.0 / 1000
    eta_residual = r_grid/8
    eta_gap = m*h/16
    eta = eta_residual + eta_gap
    # 代码检查：比较采样导数与解析构造的全域界。
    v = model(xr)
    v3 = dx(dx(dx(v, xr), xr), xr)
    assert v3.detach().abs().max().item() <= m + 1e-10
    assert e_grid <= eta + 1e-10
    assert bc == 0.0
    return dict(step=step, r_grid=r_grid, error_grid=e_grid, M=m,
                eta_residual=eta_residual, eta_gap=eta_gap, eta=eta,
                ratio_grid=eta/e_grid if e_grid else float("inf"))


x_train_base = torch.linspace(0, 1, 130)[1:-1].reshape(-1, 1)
optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
loss_history = []
records = [evaluate(model, 0)]
start = time.perf_counter()

for step in range(1, 2001):
    x_train = x_train_base.detach().clone().requires_grad_(True)
    optimizer.zero_grad(set_to_none=True)
    r = residual(model, x_train)
    loss = (r**2).mean()
    loss.backward()
    optimizer.step()
    loss_history.append(loss.item())
    if step in (100, 500, 2000):
        records.append(evaluate(model, step))

elapsed = time.perf_counter() - start
print(f"Training plus checkpoint evaluation: {elapsed:.2f} seconds")
print("step    reference-grid error    grid max residual     bound eta         M")
for row in records:
    print(f'{row["step"]:4d}    {row["error_grid"]:20.6e} '
          f'{row["r_grid"]:20.6e} {row["eta"]:14.6e} {row["M"]:10.4f}')
last = records[-1]
print("Final bound decomposition:")
print(f'  residual term R_grid/8: {last["eta_residual"]:.6e}')
print(f'  gap term M*h/16:        {last["eta_gap"]:.6e}')
print(f'  bound / grid error:    {last["ratio_grid"]:.3f}')
assert last["error_grid"] < records[0]["error_grid"]


xp = torch.linspace(0, 1, 1001).reshape(-1, 1).requires_grad_(True)
vp = model(xp)
rp = residual(model, xp)
xx = xp.detach().numpy().flatten()
vv = vp.detach().numpy().flatten()
uu = xx*(1-xx)
rr = rp.detach().numpy().flatten()

fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
axes[0, 0].plot(xx, uu, label="Exact solution")
axes[0, 0].plot(xx, vv, "--", label="PINN")
axes[0, 0].set(title="Solution after 2000 updates", xlabel="x", ylabel="u")
axes[0, 0].legend()
axes[0, 1].plot(xx, vv-uu)
axes[0, 1].set(title="Solution error on evaluation grid", xlabel="x", ylabel="v - u_exact")
axes[1, 0].plot(xx, rr)
axes[1, 0].set(title="PDE residual on evaluation grid", xlabel="x", ylabel="-v'' - 2")
axes[1, 1].semilogy(np.arange(1, 2001), loss_history)
axes[1, 1].set(title="Training residual MSE (before each update)", xlabel="Update", ylabel="Loss")
for ax in axes.flat:
    ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(OUT_DIR / "figure_01.png", dpi=150, bbox_inches="tight")
plt.close(fig)


print("intervals   grid residual term    gap correction       total bound")
for intervals in (100, 1000, 10000):
    xg = torch.linspace(0, 1, intervals+1).reshape(-1, 1).requires_grad_(True)
    rg = residual(model, xg).detach().abs().max().item()
    term_r = rg/8
    term_gap = last["M"]/(16*intervals)
    print(f"{intervals:9d} {term_r:20.6e} {term_gap:18.6e} {term_r+term_gap:18.6e}")

(OUT_DIR / "training_records.json").write_text(json.dumps(records, indent=2), encoding="utf-8")
(OUT_DIR / "final_weights.json").write_text(json.dumps({k: v.detach().tolist() for k, v in model.state_dict().items()}, indent=2), encoding="utf-8")

print("Saved results to:", OUT_DIR)
