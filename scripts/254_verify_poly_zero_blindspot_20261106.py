"""Post-exit integrity and exact-arithmetic checks of the blind-spot experiment."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'results/2026-11-06'
INPUT = BASE / 'poly_zero_derivative_blindspot.json'
OUTPUT = BASE / 'poly_zero_derivative_blindspot_readback.json'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    result = json.loads(INPUT.read_text())
    config = json.loads((BASE / 'poly_split_diagnostic/locked_config.json').read_text())
    for key, path in {
        'source_records_sha256': BASE / 'poly_split_diagnostic/records.json',
        'oracle_sha256': ROOT / 'scripts/107_polynomial_derivative_crosscheck_20261101.py',
        'split_sha256': ROOT / 'scripts/225_derivative_split_20261105.py',
        'poly_sha256': ROOT / 'scripts/74_polynomial_forcing_prototype.py',
        'script_sha256': ROOT / 'scripts/253_poly_zero_derivative_blindspot_20261106.py',
    }.items():
        if sha(path) != result[key]:
            raise ValueError(key)
    inputs = json.loads((ROOT / config['fixtures_config']).read_text())
    networks = {}
    for fixture in inputs['fixtures']:
        path = ROOT / fixture['path']
        if sha(path) != fixture['sha256']:
            raise ValueError('fixture digest mismatch')
        for net in json.loads(path.read_text())['networks']:
            if net['source_name'] in ('constant_2', 'linear_6x'):
                if net['label'] in networks:
                    raise ValueError('duplicate network')
                networks[net['label']] = (net, fixture['sha256'])
    normal = {(r['network'], tuple(r['segment']), r['order']): r
              for r in json.loads((BASE / 'poly_split_diagnostic/records.json').read_text())
              if r['pieces'] == 4 and r['location'] == 'middle'}
    keys = set()
    counts = {0: 0, 3: 0, 4: 0}
    for row in result['rows']:
        key = (row['network'], tuple(row['segment']), row['order'])
        if key in keys or row['network'] not in networks or list(key[1]) not in config['segments'] or key[2] not in counts:
            raise ValueError('unexpected or duplicate case')
        keys.add(key)
        net, digest = networks[key[0]]
        a, b = map(F, row['segment'])
        x = (a + b) / 2
        coeff = list(map(F, net['source_coefficients']))
        source = sum((v * x**j for j, v in enumerate(coeff)), F(0))
        if row['fixture_sha256'] != digest or row['source'] != net['source_name'] or F(row['x']) != x or F(row['source_value']) != source or source == 0:
            raise ValueError('fixture or source mismatch')
        # Verify the source derivative is exactly zero in the high-order checks.
        if key[2] and any(coeff[j] != 0 for j in range(key[2], len(coeff))):
            raise ValueError('higher derivative source is not zero')
        lo, hi, r80, r110 = map(F, (row['lower'], row['upper'], row['reference80'], row['reference110']))
        stable = abs(r80-r110) <= F(config['relative_stability']) * max(F(1), abs(r110))
        detected = not lo <= r110 <= hi
        expected = key[2] == 0
        if lo > hi or not stable or not row['stable'] or detected != expected or row['detected'] != detected or row['expected_detected'] != expected:
            raise ValueError('arithmetic or detection mismatch')
        if key[2] and (row['reference80'] != normal[key]['reference80'] or row['reference110'] != normal[key]['reference110']):
            raise ValueError('reference mismatch')
        counts[key[2]] += 1
    if len(networks) != 16 or counts != {0: 48, 3: 48, 4: 48} or len(keys) != 144:
        raise ValueError('incomplete run')
    if any(result[k] != v for k, v in {'networks': 16, 'expected': 144, 'records': 144, 'stable': 144, 'residual_detected': 48, 'high_order_detected': 0, 'expected_matching': 144}.items()):
        raise ValueError('summary mismatch')
    out = {'records': len(keys), 'networks': len(networks), 'residual_cases_detected': counts[0], 'higher_derivative_cases_not_detected': counts[3] + counts[4], 'records_sha256': sha(INPUT), 'scope': 'post-exit digest, completeness and exact stored-value checks; numerical point references are not a continuous proof'}
    with OUTPUT.open('x') as f:
        json.dump(out, f, indent=2)
        f.write('\n')
    print(json.dumps(out))

if __name__ == '__main__':
    main()
