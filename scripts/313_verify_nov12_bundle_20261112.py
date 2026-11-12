"""Read-only verification of the selected November 12 evidence package.

Hash the publication manifest, then rerun original readback scripts in an
isolated temporary copy. Does not rerun numerical experiments or prove bounds.
"""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "results/2026-11-12/publication_manifest.json"
CHECKS = [
    (295, "residual_gap_audit"),
    (297, "green_quadrature_audit"),
    (299, "nonnegative_residual_exact_controls"),
    (301, "sine_green_quadrature_audit"),
    (303, "polynomial_peak_audit"),
    (305, "extra_polynomial_green_audit"),
    (306, "legacy_mpmath_interval_audit"),
    (308, "mpmath_global_residual_audit"),
    (310, "mpmath_global_refinement"),
    (312, "mpmath_point_error_checks"),
]


def main():
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    files = manifest["files"]
    if not files:
        raise ValueError("empty evidence manifest")
    for relative, expected in files.items():
        path = (ROOT / relative).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file():
            raise ValueError("missing or invalid path: " + relative)
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError("hash mismatch: " + relative)
    checks = []
    with tempfile.TemporaryDirectory(prefix="pde-nov12-readback-") as temp:
        target = Path(temp)
        for relative in files:
            dest = target / relative
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes((ROOT / relative).read_bytes())
        env = dict(os.environ)
        env.pop("PYTHONOPTIMIZE", None)  # Original checkers include assert statements.
        env["PYTHONDONTWRITEBYTECODE"] = "1"
        for number, directory in CHECKS:
            scripts = list(target.glob(f"scripts/{number}_*20261112.py"))
            if len(scripts) != 1:
                raise ValueError(f"expected one checker for {number}")
            output = target / "results/2026-11-12" / directory / "independent_readback.json"
            expected = output.read_bytes()
            output.unlink()  # Only the disposable copy; archived evidence is untouched.
            run = subprocess.run(
                [sys.executable, str(scripts[0])], cwd=target, env=env,
                capture_output=True, text=True, timeout=120,
            )
            if run.returncode:
                raise RuntimeError(f"checker {number} failed:\n{run.stdout}\n{run.stderr}")
            if not output.is_file() or output.read_bytes() != expected:
                raise ValueError(f"checker {number} readback differs from archive")
            checks.append({"script": scripts[0].name, "readback_equal": True})
    print(json.dumps({
        "manifest_files_checked": len(files), "readbacks": checks,
        "python": sys.version.split()[0],
        "scope": "file hashes and ten original stored-result readbacks; no numerical-kernel rerun, Windows test, or continuous-domain proof",
    }, indent=2))


if __name__ == "__main__":
    main()
