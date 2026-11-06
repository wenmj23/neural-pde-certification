"""Replay selected 251–254 experiments in a fresh, non-overwriting directory."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir', required=True, type=Path, help='new, nonexistent output directory')
    args = parser.parse_args()
    work = args.workdir.resolve()
    if work.exists():
        raise FileExistsError(work)
    # Reject an ancestor of the repository; never modify repository files.
    if work == ROOT or work in ROOT.parents:
        raise ValueError('workdir cannot be the repository or its ancestor')
    work.mkdir(parents=True)
    try:
        # Include Python modules and input configurations as a private workspace.
        shutil.copytree(ROOT/'scripts', work/'scripts', ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
        shutil.copytree(ROOT/'configs', work/'configs', ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
        for name in ('locked_config.json', 'summary.json'):
            rel = Path('results/2026-11-05/independent_poly_split')/name
            dst = work/rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT/rel, dst)
        for script in ('251_poly_split_diagnostic_20261106.py','252_verify_poly_split_diagnostic_20261106.py','253_poly_zero_derivative_blindspot_20261106.py','254_verify_poly_zero_blindspot_20261106.py'):
            subprocess.run([sys.executable, str(work/'scripts'/script)], cwd=work, check=True)
        paths = ('results/2026-11-06/poly_split_diagnostic/records.json','results/2026-11-06/poly_split_diagnostic/faults.json','results/2026-11-06/poly_zero_derivative_blindspot.json')
        comparisons = {p: {'saved_sha256': sha(ROOT/p), 'replay_sha256': sha(work/p), 'identical': sha(ROOT/p)==sha(work/p)} for p in paths}
        result = {'status': 'checks_passed', 'bytewise_comparison': comparisons, 'scope': 'selected numerical point diagnostics only; not an independent mathematical proof; different platforms may produce different bytes'}
        (work/'replay_report.json').write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps(result, indent=2), flush=True)
    except Exception:
        # Preserve failed work for diagnosis; no existing results are replaced.
        print(f'Replay incomplete; files kept at {work}', file=sys.stderr)
        raise

if __name__ == '__main__':
    main()
