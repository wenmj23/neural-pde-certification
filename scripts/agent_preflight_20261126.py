"""Portable, non-destructive readiness check and optional public CLI smoke run."""
import argparse
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--smoke-out', type=Path, help='NEW directory for a single smoke run')
    args = parser.parse_args()
    paths = ['AGENTS.md', 'START_HERE.md', 'HANDOFF.md', 'NEXT_STEPS.md',
             'scripts/certify.py', 'configs/examples/shallow_sine.json',
             'configs/examples/deep_sine.json', 'requirements-certifier.txt']
    versions = {}
    for package in ('python-flint', 'mpmath'):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    hashes = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
              for p in paths if (ROOT / p).is_file()}
    try:
        git = subprocess.run(['git', 'status', '--short', '--branch'], cwd=ROOT,
                             capture_output=True, text=True, timeout=15)
        git_info = {'returncode': git.returncode, 'stdout': git.stdout, 'stderr': git.stderr}
    except (OSError, subprocess.TimeoutExpired) as exc:
        git_info = {'error': str(exc)}
    data = {'utc': datetime.now(timezone.utc).isoformat(), 'python': sys.version,
            'platform': platform.platform(), 'packages': versions,
            'missing_files': [p for p in paths if p not in hashes],
            'file_sha256': hashes, 'git_status': git_info,
            'scope': 'Readiness and optional existing-case smoke run; not a mathematical proof or agent scheduler.'}
    if data['missing_files'] or any(v is None for v in versions.values()):
        print(json.dumps(data, ensure_ascii=False, indent=2))
        return 1
    if args.smoke_out:
        dest = args.smoke_out.expanduser().resolve()
        if dest.exists():
            raise FileExistsError(f'Output path already exists: {dest}')
        command = [sys.executable, str(ROOT / 'scripts/certify.py'), '--case',
                   str(ROOT / 'configs/examples/shallow_sine.json'), '--out', str(dest),
                   '--degree', '2', '--mesh', 'uniform', '--goal', 'decision',
                   '--error-tolerance', '1e-3', '--max-intervals', '32']
        try:
            run = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=90)
            data['smoke'] = {'command': command, 'returncode': run.returncode,
                             'stdout': run.stdout, 'stderr': run.stderr,
                             'result_exists': (dest / 'result.json').is_file()}
        except subprocess.TimeoutExpired:
            data['smoke'] = {'command': command, 'timeout': True}
        if 'returncode' not in data['smoke'] or data['smoke']['returncode'] not in (0, 2) or not data['smoke']['result_exists']:
            print(json.dumps(data, ensure_ascii=False, indent=2))
            return 2
    print(json.dumps(data, ensure_ascii=False, indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
