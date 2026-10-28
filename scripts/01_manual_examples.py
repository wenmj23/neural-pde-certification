# -*- coding: utf-8 -*-
"""Experiment 01: extracted from the executed notebook. See README for assumptions."""
from pathlib import Path
import json
import matplotlib
matplotlib.use("Agg")
OUT_DIR = Path(__file__).resolve().parents[1] / "results" / "generated" / "01_manual_examples"
OUT_DIR.mkdir(parents=True, exist_ok=True)

import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0.0, 1.0, 20001)
u_exact = x * (1.0 - x)

# 案例 A：平缓误差。解析最大误差为 0.1/4。
# 案例 B：振荡误差。二阶导数使振幅乘以 (4*pi)^2。
cases = [
    dict(name="A: smooth arch", error=0.1*x*(1-x),
         epsilon=0.2, true_max=0.025),
    dict(name="B: oscillatory", error=0.01*np.sin(4*np.pi*x),
         epsilon=0.01*(4*np.pi)**2, true_max=0.01),
]

print("case                 residual bound    error bound    true max error    ratio")
for case in cases:
    case["eta"] = case["epsilon"] / 8
    case["ratio"] = case["eta"] / case["true_max"]
    assert case["eta"] >= case["true_max"] - 1e-14
    assert np.max(np.abs(case["error"])) <= case["eta"] + 1e-14
    assert np.allclose(case["error"][[0, -1]], 0.0, atol=1e-14)
    print(f'{case["name"]:20s} {case["epsilon"]:14.6f} '
          f'{case["eta"]:14.6f} {case["true_max"]:17.6f} {case["ratio"]:8.3f}')


fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
for ax, case in zip(axes, cases):
    q = case["epsilon"] / 2 * x * (1-x)
    ax.fill_between(x, -q, q, alpha=0.16, color="tab:blue", label="Guaranteed envelope")
    ax.plot(x, q, "--", color="tab:blue", linewidth=1)
    ax.plot(x, -q, "--", color="tab:blue", linewidth=1)
    ax.plot(x, case["error"], color="tab:red", label="Actual error e(x)")
    ax.axhline(0, color="0.6", linewidth=0.7)
    ax.set(title=case["name"], xlabel="Position x", ylabel="Error")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8)
fig.suptitle("The same valid inequality can be tight or conservative")
fig.tight_layout()
fig.savefig(OUT_DIR / "figure_01.png", dpi=150, bbox_inches="tight")
plt.close(fig)


k = 20
amplitude = 0.01
sample_x = np.linspace(0.0, 1.0, 11)
def residual(z):
    return amplitude * (k*np.pi)**2 * np.sin(k*np.pi*z)

sample_max = np.max(np.abs(residual(sample_x)))
analytic_epsilon = amplitude * (k*np.pi)**2
print(f"Residual maximum observed at 11 points: {sample_max:.3e}")
print(f"Analytic global residual maximum:       {analytic_epsilon:.6f}")
print(f"Incorrect plug-in bound from sampling:  {sample_max/8:.3e}")
print(f"Valid bound from analytic residual:     {analytic_epsilon/8:.6f}")
print(f"True maximum solution error:            {amplitude:.6f}")
assert sample_max/8 < amplitude
assert analytic_epsilon/8 >= amplitude

fig, ax = plt.subplots(figsize=(10, 3.6))
ax.plot(x, np.abs(residual(x)), label="Absolute residual between samples")
ax.scatter(sample_x, np.abs(residual(sample_x)), color="tab:red", zorder=3,
           label="11 equally spaced check points")
ax.set(xlabel="Position x", ylabel="Absolute residual",
       title="A sampled maximum can miss the continuous maximum")
ax.grid(alpha=0.2)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUT_DIR / "figure_02.png", dpi=150, bbox_inches="tight")
plt.close(fig)

print("Saved results to:", OUT_DIR)
