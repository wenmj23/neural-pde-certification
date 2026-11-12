"""Retrospective finite-grid audit of residual-bound looseness.

The NumPy forward/derivative implementation and linearly reconstructed Green
integral are independent diagnostics. They are not validated intervals.
Previously reported bounds/references are read, never recomputed or replaced.
"""
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import os
import platform
import sys
import tempfile
import time

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CFG = ROOT / "configs/residual_gap_audit_2026-11-12.json"
DEST = ROOT / "results/2026-11-12/residual_gap_audit"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def value_and_residual(net, x):
    a = x.reshape(-1, 1)
    a1 = np.ones_like(a)
    a2 = np.zeros_like(a)
    for layer in net["layers"]:
        w = np.asarray(layer["weight"], dtype=np.float64)
        b = np.asarray(layer["bias"], dtype=np.float64)
        z, z1, z2 = a @ w.T + b, a1 @ w.T, a2 @ w.T
        if layer["activation"] == "tanh":
            a = np.tanh(z)
            s = 1 - a*a
            a1 = s*z1
            a2 = s*z2 - 2*a*s*z1*z1
        elif layer["activation"] == "identity":
            a, a1, a2 = z, z1, z2
        else:
            raise ValueError("unsupported activation")
    if a.shape[1] != 1:
        raise ValueError("expected one output")
    n, n1, n2 = a[:, 0], a1[:, 0], a2[:, 0]
    h = x*(1-x)
    v = h*n
    v2 = h*n2 + 2*(1-2*x)*n1 - 2*n
    f = np.polynomial.polynomial.polyval(x, net["source_coefficients"])
    u = np.zeros_like(x)
    for j, c in enumerate(net["source_coefficients"]):
        u += c*(x-x**(j+2))/((j+1)*(j+2))
    return v-u, -v2-f


def green_linear_at_nodes(y):
    """Integrate Green's kernel against the piecewise-linear nodal data."""
    n = len(y)-1
    h = 1.0/n
    x = np.linspace(0.0, 1.0, n+1)
    area = h*(y[:-1]+y[1:])/2
    moment = x[:-1]*area + h*h*(y[:-1]+2*y[1:])/6
    left = np.concatenate(([0.0], np.cumsum(moment)))
    right_cell = area-moment
    right = np.concatenate((np.cumsum(right_cell[::-1])[::-1], [0.0]))
    return (1-x)*left+x*right


def high_precision(net, index, grid, digits):
    with mp.workdps(digits):
        x = mp.mpf(index)/grid

        def candidate(s):
            a = [s]
            for layer in net["layers"]:
                z = [sum((mp.mpf(w)*v for w, v in zip(row, a)), mp.mpf(b))
                     for row, b in zip(layer["weight"], layer["bias"])]
                if layer["activation"] == "tanh":
                    a = [mp.tanh(q) for q in z]
                elif layer["activation"] == "identity":
                    a = z
                else:
                    raise ValueError("unsupported activation")
            return s*(1-s)*a[0]

        coeffs = net["source_coefficients"]
        u = sum((mp.mpf(c)*(x-x**(j+2))/((j+1)*(j+2))
                 for j, c in enumerate(coeffs)), mp.mpf(0))
        f = sum((mp.mpf(c)*x**j for j, c in enumerate(coeffs)), mp.mpf(0))
        return mp.nstr(candidate(x)-u, digits-5), mp.nstr(-mp.diff(candidate, x, 2)-f, digits-5)


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    config = json.loads(CFG.read_text())
    fixture, bp = ROOT/config["fixture_path"], ROOT/config["baseline_path"]
    if sha(fixture) != config["fixture_sha256"] or sha(bp) != config["baseline_sha256"]:
        raise ValueError("trusted input hash mismatch")
    networks = json.loads(fixture.read_text())["networks"]
    old = json.loads(bp.read_text())
    if len(networks) != 12 or len({n["label"] for n in networks}) != 12:
        raise ValueError("network manifest mismatch")
    baselines = {(r["network"], r["method"]): r for r in old
                 if r["intervals"] == config["reference_mesh"]}
    methods = ["global", "spatial", "signed_quadratic", "signed_cubic"]
    if len(baselines) != 48 or any((n["label"], m) not in baselines for n in networks for m in methods):
        raise ValueError("baseline manifest mismatch")
    rows, spots, synthetic = [], [], []
    t0 = time.monotonic()
    for net in networks:
        label = net["label"]
        saved = {m: baselines[(label, m)] for m in methods}
        midpoint = float((Fraction(saved["global"]["reference_lower"]) +
                          Fraction(saved["global"]["reference_upper"]))/2)
        bounds = {m: float(Fraction(saved[m]["upper"])) for m in methods}
        if midpoint <= 0:
            raise ValueError("reference midpoint must be positive")
        for grid in config["grids"]:
            x = np.linspace(0, 1, grid+1)
            e, r = value_and_residual(net, x)
            zg, ag = green_linear_at_nodes(r), green_linear_at_nodes(np.abs(r))
            err = float(np.max(np.abs(e)))
            absprop = float(np.max(ag))
            if not all(np.isfinite(q) for q in [err, absprop, *bounds.values()]):
                raise ArithmeticError("nonfinite diagnostic")
            rows.append({"network": label, "grid_intervals": grid,
                         "grid_max_abs_error": err,
                         "grid_max_abs_residual": float(np.max(np.abs(r))),
                         "grid_min_residual": float(np.min(r)),
                         "grid_max_residual": float(np.max(r)),
                         "both_signs_on_grid": bool(np.min(r)<0<np.max(r)),
                         "green_signed_grid_max": float(np.max(np.abs(zg))),
                         "green_abs_grid_max": absprop,
                         "green_direct_max_difference": float(np.max(np.abs(zg-e))),
                         "sampled_global_over_sampled_error": float(np.max(np.abs(r))/8/err),
                         "sampled_abs_green_over_sampled_error": absprop/err,
                         "sampled_global_over_sampled_abs_green": float(np.max(np.abs(r))/8/absprop),
                         "reference_midpoint": midpoint, "saved_64_bounds": bounds,
                         "saved_global_over_reference": bounds["global"]/midpoint,
                         "saved_spatial_over_reference": bounds["spatial"]/midpoint,
                         "sample_error_below_saved_cubic": bool(err <= bounds["signed_cubic"]),
                         "sample_abs_green_below_saved_spatial": bool(absprop <= bounds["spatial"])})
        grid = max(config["grids"])
        x = np.linspace(0, 1, grid+1)
        e, r = value_and_residual(net, x)
        indices = {"max_abs_error": int(np.argmax(np.abs(e))),
                   "max_residual": int(np.argmax(r)),
                   "min_residual": int(np.argmin(r)),
                   "one_third": grid//3}
        for kind in config["high_precision_spots"]:
            idx = indices[kind]
            v80 = high_precision(net, idx, grid, config["high_precision_digits"][0])
            v110 = high_precision(net, idx, grid, config["high_precision_digits"][1])
            with mp.workdps(130):
                sdiff = [mp.nstr(abs(mp.mpf(a)-mp.mpf(b)), 8) for a, b in zip(v80,v110)]
            spots.append({"network":label,"kind":kind,"grid_index":idx,"grid_intervals":grid,
                          "error_80":v80[0],"error_110":v110[0],
                          "residual_80":v80[1],"residual_110":v110[1],
                          "precision_abs_differences":sdiff,
                          "float_abs_error_difference":abs(float(e[idx])-float(v110[0])),
                          "float_residual_difference":abs(float(r[idx])-float(v110[1]))})
        print(label, "complete", flush=True)
    grid = max(config["grids"])
    x = np.linspace(0, 1, grid+1)
    for k in config["manufactured_frequencies"]:
        e = np.sin(k*np.pi*x)
        r = (k*np.pi)**2*e
        z = green_linear_at_nodes(r)
        a = green_linear_at_nodes(np.abs(r))
        synthetic.append({"frequency":k,"grid_intervals":grid,
                          "exact_error_max":1.0,"exact_global_ratio":float((k*np.pi)**2/8),
                          "green_signed_max":float(np.max(np.abs(z))),
                          "sampled_abs_green_max":float(np.max(a)),
                          "max_point_difference":float(np.max(np.abs(z-e)))})
    if len(rows) != 60 or len(spots) != 48 or len(synthetic) != 9:
        raise ValueError("output manifest mismatch")
    base = DEST.parent
    base.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".residual_gap_stage_", dir=base))
    (stage/"locked_config.json").write_bytes(CFG.read_bytes())
    for name, data in [("records.json", rows), ("spot_checks.json", spots), ("synthetic.json", synthetic)]:
        (stage/name).write_text(json.dumps(data, indent=2, allow_nan=False)+"\n")
    summaries = {"expected_records":60,"records":len(rows),"expected_spots":48,"spots":len(spots),
                 "expected_synthetic":9,"synthetic":len(synthetic),
                 "fixture_sha256":sha(fixture),"baseline_sha256":sha(bp),
                 "script_sha256":sha(Path(__file__)),"config_sha256":sha(CFG),
                 "result_sha256":{name:sha(stage/name) for name in ["records.json","spot_checks.json","synthetic.json"]},
                 "elapsed_seconds":time.monotonic()-t0,"python":sys.version,
                 "platform":platform.platform(),"numpy_version":np.__version__,"mpmath_version":mp.__version__,
                 "scope":config["scope"]}
    (stage/"summary.json").write_text(json.dumps(summaries, indent=2, allow_nan=False)+"\n")
    os.rename(stage, DEST)
    print(json.dumps(summaries), flush=True)


if __name__ == "__main__":
    main()
