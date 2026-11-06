"""Read back 255 output: hashes, full case manifest and exact arithmetic."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import factorial
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'results/2026-11-06/multiorder_fault_diagnostic'
OUT = ROOT / 'results/2026-11-06/multiorder_fault_diagnostic_readback.json'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    if OUT.exists():
        raise FileExistsError(OUT)
    s = json.loads((DATA / 'summary.json').read_text())
    c = json.loads((DATA / 'locked_config.json').read_text())
    for key, path in {'records_sha256': DATA/'records.json', 'blindspot_sha256': DATA/'blindspot.json', 'config_sha256': DATA/'locked_config.json', 'fixture_config_sha256': ROOT/c['fixtures_config'], 'poly_sha256': ROOT/'scripts/74_polynomial_forcing_prototype.py', 'oracle_sha256': ROOT/'scripts/107_polynomial_derivative_crosscheck_20261101.py', 'script_sha256': ROOT/'scripts/255_multiorder_fault_diagnostic_20261106.py'}.items():
        if sha(path) != s[key]:
            raise ValueError('digest mismatch: ' + key)
    fixtures = json.loads((ROOT/c['fixtures_config']).read_text())
    networks = {}
    for fixture in fixtures['fixtures']:
        path = ROOT / fixture['path']
        if sha(path) != fixture['sha256']:
            raise ValueError('fixture hash mismatch')
        for net in json.loads(path.read_text())['networks']:
            if net['label'] in networks:
                raise ValueError('duplicate network')
            networks[net['label']] = (net, fixture['sha256'])
    points = [F(x) for x in c['points']]
    kinds = ['normal'] + c['faults']
    expected = set(product(networks, points, c['orders'], kinds))
    seen = set()
    observations = {}
    counts = Counter()
    group = Counter()
    for row in json.loads((DATA/'records.json').read_text()):
        key = row['network'], F(row['x']), row['order'], row['kind']
        if key not in expected or key in seen:
            raise ValueError('invalid or duplicate case')
        seen.add(key)
        net, digest = networks[key[0]]
        if row['fixture_sha256'] != digest or row['source'] != net['source_name']:
            raise ValueError('fixture mismatch')
        src = sum((F(v) * factorial(i) / factorial(i-key[2]) * key[1]**(i-key[2]) for i,v in enumerate(net['source_coefficients']) if i >= key[2]), F(0))
        if F(row['source_derivative']) != src:
            raise ValueError('source derivative mismatch')
        lo, hi, a, b = map(F, (row['lower'], row['upper'], row['reference80'], row['reference110']))
        stable = abs(a-b) <= F(c['relative_stability']) * max(F(1), abs(b))
        contained = lo <= b <= hi
        if lo > hi or not stable or not row['stable'] or row['contained'] != contained:
            raise ValueError('stored arithmetic mismatch')
        if contained != (key[3] == 'normal' or src == 0):
            raise ValueError('unexpected injected fault result')
        observations[key] = row
        counts[(key[3], 'contained' if contained else 'detected')] += 1
        group[(row['source'], key[2], key[3], 'zero' if src == 0 else 'nonzero', 'contained' if contained else 'detected')] += 1
    if seen != expected or len(networks) != 24:
        raise ValueError('incomplete case manifest')
    for name, x, order in product(networks, points, c['orders']):
        normal = observations[(name,x,order,'normal')]
        for kind in c['faults']:
            other = observations[(name,x,order,kind)]
            if (normal['reference80'], normal['reference110']) != (other['reference80'], other['reference110']):
                raise ValueError('reference reuse mismatch')
    blind = json.loads((DATA/'blindspot.json').read_text())
    extra = F(c['extra_point'])
    expected_blind = set(product(networks, points+[extra]))
    blind_seen = set()
    bcounts = Counter()
    for row in blind:
        key = row['network'], F(row['x'])
        if key in blind_seen or key not in expected_blind or row['fixture_sha256'] != networks[key[0]][1]:
            raise ValueError('invalid blindspot case')
        blind_seen.add(key)
        delta = F(1)
        for p in points:
            delta *= key[1]-p
        lo,hi,a,b = map(F,(row['lower'],row['upper'],row['reference80'],row['reference110']))
        contained = lo <= b <= hi
        if F(row['delta']) != delta or row['contained'] != contained or contained != (key[1] in points) or abs(a-b) > F(c['relative_stability'])*max(F(1),abs(b)):
            raise ValueError('blindspot arithmetic')
        if key[1] in points:
            normal = observations[(key[0],key[1],0,'normal')]
            if (row['lower'],row['upper'],row['reference80'],row['reference110']) != (normal['lower'],normal['upper'],normal['reference80'],normal['reference110']):
                raise ValueError('sampled blindspot differs from normal')
        bcounts['sample_contained' if key[1] in points else 'extra_detected'] += 1
    if blind_seen != expected_blind:
        raise ValueError('incomplete blindspot cases')
    totals = {'records': len(seen), 'blindspot_records': len(blind_seen), 'stable': len(seen), 'normal_contained': counts['normal','contained'], 'nonzero_fault_detected': sum(counts[k,'detected'] for k in c['faults']), 'zero_fault_contained': sum(counts[k,'contained'] for k in c['faults']), 'sampled_blindspot_contained': bcounts['sample_contained'], 'extra_blindspot_detected': bcounts['extra_detected']}
    if any(s[k] != v for k,v in totals.items()):
        raise ValueError('summary mismatch')
    output = dict(totals, records_sha256=sha(DATA/'records.json'), blindspot_sha256=sha(DATA/'blindspot.json'), groups=[{'source':k[0], 'order':k[1], 'kind':k[2], 'source_derivative':k[3], 'outcome':k[4], 'count':v} for k,v in sorted(group.items())], scope='readback of stored numeric diagnostics and exact arithmetic; not a continuous proof')
    with OUT.open('x') as f:
        json.dump(output,f,indent=2)
        f.write('\n')
    print(json.dumps(totals))

if __name__ == '__main__':
    main()
