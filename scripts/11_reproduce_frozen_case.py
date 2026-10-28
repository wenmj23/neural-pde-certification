"""Reload the canonical fixed PINN and reproduce reference/certificate outputs.

No optimizer or training is run. Ordinary float64 baseline only; the exact
solution is used exclusively for reference evaluation and implementation checks.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np
import torch
from torch import nn
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from certification_core import NeuralResidual, certificate

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CASE = ROOT / 'cases/poisson_1d_sine_k1_seed7'


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FrozenPINN(nn.Module):
    def __init__(self, weights):
        super().__init__()
        self.hidden = nn.Linear(1, 16, dtype=torch.float64)
        self.output = nn.Linear(16, 1, dtype=torch.float64)
        self.load_state_dict({k: torch.tensor(v, dtype=torch.float64)
                              for k, v in weights.items()}, strict=True)
        self.eval()

    def forward(self, x):
        return x * (1 - x) * self.output(torch.tanh(self.hidden(x)))


def load_case(directory):
    case = json.loads((directory / 'case.json').read_text())
    if case['schema_version'] != 1:
        raise ValueError('Unsupported case schema')
    p, n = case['problem'], case['network']
    if (p['operator'] != '-u_dd' or p['domain'] != [0., 1.]
            or p['boundary'] != {'type': 'Dirichlet', 'values': [0., 0.]}
            or p['frequency'] != 1 or p['source'] != 'pi^2*sin(pi*x)'
            or p['exact_solution_for_reference_only'] != 'sin(pi*x)'
            or n['architecture'] != [1, 16, 1] or n['activation'] != 'tanh'
            or n['output_transform'] != 'x*(1-x)*raw_network(x)'
            or n['parameter_dtype'] != 'float64'):
        raise ValueError('This reproduction script only supports the stated frozen case')
    weights_path = directory / n['weights_file']
    provenance_path = directory / case['provenance_file']
    if sha256(weights_path) != n['weights_sha256']:
        raise ValueError('Frozen weights have changed')
    if sha256(provenance_path) != case['provenance_sha256']:
        raise ValueError('Frozen provenance has changed')
    if 'baseline_file' in case:
        if sha256(directory / case['baseline_file']) != case['baseline_sha256']:
            raise ValueError('Frozen baseline has changed')
    weights = json.loads(weights_path.read_text())
    model = FrozenPINN(weights)
    assert sum(p.numel() for p in model.parameters()) == n['parameter_count'] == 49
    assert all(torch.isfinite(p).all() for p in model.parameters())
    params = (np.array(weights['hidden.weight']).flatten(),
              np.array(weights['hidden.bias']),
              np.array(weights['output.weight']).flatten(), weights['output.bias'][0])
    residual = NeuralResidual(params, p['frequency'], case['case_id'])
    return case, model, residual


def solution_numpy(x, params):
    w, b, a, c = params
    return x * (1 - x) * (c + (np.tanh(x[..., None] * w + b) * a).sum(axis=-1))


def cross_check(model, problem):
    x = torch.linspace(0, 1, 257, dtype=torch.float64).reshape(-1, 1).requires_grad_(True)
    v = model(x)
    def dx(y):
        return torch.autograd.grad(y, x, torch.ones_like(y), create_graph=True)[0]
    r = -dx(dx(v)) - np.pi ** 2 * torch.sin(np.pi * x)
    xnp = x.detach().numpy().flatten()
    v_diff = np.max(abs(v.detach().numpy().flatten() - solution_numpy(xnp, problem.params)))
    r_diff = np.max(abs(r.detach().numpy().flatten() - problem.residual(xnp)))
    assert v_diff < 2e-13
    assert r_diff < 2e-11
    boundary = model(torch.tensor([[0.], [1.]], dtype=torch.float64)).detach().numpy().flatten()
    assert np.array_equal(boundary, np.zeros(2))
    return dict(points=257, max_solution_difference=float(v_diff),
                max_residual_difference=float(r_diff), boundary_values=boundary.tolist(),
                scope='Independent evaluation and autograd cross-check; not a rounding-error proof')


def compare_number(actual, expected):
    if not np.isclose(actual, expected, rtol=5e-9, atol=5e-14):
        raise AssertionError(f'Frozen baseline mismatch: {actual} versus {expected}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', type=Path, default=DEFAULT_CASE)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    torch.set_num_threads(1)
    case_dir = args.case.resolve()
    case, model, problem = load_case(case_dir)
    weights_digest_before = sha256(case_dir / case['network']['weights_file'])
    cross = cross_check(model, problem)
    reference = case['reference']
    x = np.linspace(0, 1, reference['points'])
    v = solution_numpy(x, problem.params)
    exact = np.sin(np.pi * x)
    error = v - exact
    peak = int(np.argmax(abs(error)))
    grid_error = float(abs(error[peak]))
    compare_number(grid_error, reference['error_grid'])
    rows = []
    for intervals in case['reproduction']['interval_budgets']:
        nodes = np.linspace(0, 1, intervals + 1)
        for method in ('constant', 'linear'):
            # This call does not receive the exact solution or reference error.
            r = certificate(problem, nodes, method)
            assert r['upper'] >= grid_error - 5e-14
            rows.append(dict(intervals=intervals, method=method, lower=r['lower'],
                             upper=r['upper'], reconstructed_maximum=r['Z'], remainder=r['D'],
                             ratio_to_grid_error=r['upper'] / grid_error))
    for old in reference['historical_certificates_1000_intervals']:
        new = next(r for r in rows if r['intervals'] == 1000 and r['method'] == old['method'])
        for field in ('lower', 'upper', 'remainder'):
            compare_number(new[field], old[field])
    baseline_path = case_dir / case.get('baseline_file', 'baseline.json')
    baseline_checked = baseline_path.exists()
    if baseline_checked:
        old = json.loads(baseline_path.read_text())
        compare_number(grid_error, old['reference']['error_grid'])
        for new, previous in zip(rows, old['certificates'], strict=True):
            assert (new['intervals'], new['method']) == (previous['intervals'], previous['method'])
            for field in ('lower', 'upper', 'reconstructed_maximum', 'remainder'):
                compare_number(new[field], previous[field])
    assert sha256(case_dir / case['network']['weights_file']) == weights_digest_before
    out = args.out or ROOT / 'results/generated/11_frozen_main_case' / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True, exist_ok=False)
    result = dict(case_id=case['case_id'], weights_sha256=weights_digest_before,
                  parameter_count=49, retrained=False,
                  reference=dict(kind='observed_dense_grid_maximum', points=len(x),
                                 error_grid=grid_error, max_error_grid_location=float(x[peak]),
                                 exact_solution_used_by_certificate=False,
                                 caveat='Grid maximum is not the exact continuous supremum'),
                  certificates=rows, torch_numpy_cross_check=cross,
                  frozen_baseline_checked=baseline_checked,
                  arithmetic='ordinary float64; machine-level rounding audit pending')
    (out / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    environment = dict(python=sys.version, platform=platform.platform(), numpy=np.__version__,
                       torch=torch.__version__, matplotlib=matplotlib.__version__,
                       source_hashes={str(p.relative_to(ROOT)): sha256(p) for p in
                                      [Path(__file__), ROOT / 'scripts/certification_core.py']},
                       case_file_sha256=sha256(case_dir / 'case.json'),
                       provenance_sha256=sha256(case_dir / case['provenance_file']))
    (out / 'environment.json').write_text(json.dumps(environment, indent=2) + '\n')
    residual_values = problem.residual(x)
    with (out / 'reference_curve.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['x', 'pinn_value', 'reference_exact_solution', 'signed_error', 'residual'])
        writer.writerows(zip(x, v, exact, error, residual_values))
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    axes[0].plot(x, exact, label='Exact solution: sin(pi x)', color='#277da1', lw=2)
    axes[0].plot(x[::60], v[::60], '.', label='Frozen PINN', color='#f3722c', ms=3)
    axes[0].set(title='Frozen main case: 1D Poisson, k=1, seed=7', ylabel='Solution')
    axes[0].legend()
    axes[1].plot(x, abs(error), color='#43aa8b', label='Absolute error on reference grid')
    axes[1].plot(x[peak], abs(error[peak]), 'o', color='#f3722c', ms=5)
    axes[1].set(xlabel='x', ylabel='Absolute error', title=f'Observed grid maximum: {grid_error:.7e}')
    axes[1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    axes[1].legend()
    for ax in axes: ax.grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(out / 'frozen_case.png', dpi=160)
    plt.close(fig)
    print(json.dumps(result, ensure_ascii=False, indent=2), flush=True)
    print('RESULT_DIRECTORY', out, flush=True)


if __name__ == '__main__': main()
