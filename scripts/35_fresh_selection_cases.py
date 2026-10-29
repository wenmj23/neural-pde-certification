"""Freeze six unselected seed-37 networks after selector code is locked.

Training provides test inputs only. Independent continuous references are
computed afterwards, outside all selector timing and decision paths.
"""
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import time

import torch

ROOT = Path(__file__).resolve().parents[1]
CANONICAL = ROOT/'configs/frozen_networks_selection_fresh_2026-10-29.json'


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def main():
    if CANONICAL.exists():
        raise FileExistsError('Preserve existing fresh fixture')
    training = load('training22', ROOT/'scripts/22_freeze_depth_networks.py')
    reference = load('reference24', ROOT/'scripts/24_depth_continuous_reference.py')
    out = ROOT/'results/generated/35_fresh_selection_cases'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True, exist_ok=False)
    names = ['scripts/arb_degree_selection.py', 'scripts/33_degree_selection_study.py']
    locked = {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names}
    sources = [Path(__file__), *[ROOT/n for n in names],
               *[ROOT/'scripts'/n for n in ('22_freeze_depth_networks.py', '24_depth_continuous_reference.py',
                 '18_continuous_reference_crosscheck.py', 'interval_certification.py')]]
    snapshot = out/'source_snapshot'; snapshot.mkdir()
    for p in sources: (snapshot/p.name).write_bytes(p.read_bytes())
    config = dict(status='running', seed=37, depths=[1, 2, 3], frequencies=[1, 5], width=16,
                  steps=2000, optimizer='Adam', learning_rate=.01, collocation_points=128,
                  dtype='float64', threads=1, torch=torch.__version__, python=platform.python_version(),
                  locked_selection_sources=locked, reference_run=str(out.relative_to(ROOT)),
                  selection='All six final outputs retained. No observed/true error field used to select or remove cases.',
                  source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    def save(name, value):
        p = out/name; p.write_text(json.dumps(value, indent=2)+'\n')
        if json.loads(p.read_text()) != value: raise OSError('Readback mismatch')
    save('config.json', config)
    torch.set_default_dtype(torch.float64); torch.set_num_threads(1)
    rows = []
    try:
        for depth in (1, 2, 3):
            for k in (1, 5):
                torch.manual_seed(37); model = training.PINN(depth)
                optimizer = torch.optim.Adam(model.parameters(), lr=.01)
                base = torch.linspace(0, 1, 130)[1:-1].reshape(-1, 1)
                history = []; started = time.perf_counter()
                for step in range(1, 2001):
                    x = base.detach().clone().requires_grad_(True)
                    optimizer.zero_grad(set_to_none=True)
                    loss = training.residual(model, x, k).square().mean()
                    if not torch.isfinite(loss): raise ArithmeticError('Nonfinite training loss')
                    loss.backward(); optimizer.step()
                    if step in (1, 100, 500, 1000, 2000):
                        history.append(dict(step=step, loss_before_update=float(loss.detach())))
                layers = [dict(weight=l.weight.detach().tolist(), bias=l.bias.detach().tolist(),
                               activation='tanh' if i < depth else 'identity') for i, l in enumerate(model.layers)]
                row = dict(label=f'depth{depth}_k{k}_seed37', depth=depth, width=16, frequency=k,
                           seed=37, parameter_count=sum(p.numel() for p in model.parameters()),
                           layers=layers, training_history=history, training_wall_s=time.perf_counter()-started)
                rows.append(row); save('networks.json', dict(config=config, networks=rows))
                print('NEW_FROZEN', row['label'], 'seconds', round(row['training_wall_s'], 2), flush=True)
        fixture_config = {k:v for k, v in config.items() if k not in ('status', 'source_sha256')}
        fixture = dict(config=fixture_config, networks=rows)
        save('networks.json', fixture)
        CANONICAL.write_bytes((out/'networks.json').read_bytes())
        config['fixture_sha256'] = hashlib.sha256(CANONICAL.read_bytes()).hexdigest()
        refs = {}
        for row in rows:
            ref = reference.MPDeepReference(row).maximum(relative=F(1, 10**8), absolute=F(1, 10**14), cap=5000)
            refs[row['label']] = ref; save('references.json', refs)
            print('NEW_REFERENCE', row['label'], 'converged', ref['converged'], 'seconds', round(ref['elapsed_s'], 2), flush=True)
            if not ref['converged']: raise ArithmeticError('Reference budget exhausted')
        for name, sha in locked.items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != sha:
                raise RuntimeError('Selection source changed during fresh-case generation')
        config.update(status='complete', frozen_cases=len(rows), continuous_references=len(refs))
        save('config.json', config)
    except Exception as error:
        config.update(status='failed', failure=repr(error)); save('config.json', config); raise
    print('COMPLETED', out, 'six unselected fresh cases and references', flush=True)


if __name__ == '__main__': main()
