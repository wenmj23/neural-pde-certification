"""Readback the selected point-recheck records against original segment cases."""
from fractions import Fraction as F
from collections import Counter
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/2026-11-06';SRC=BASE/'low_order_interval_diagnostic/records.json';DATA=BASE/'low_order_segment_masking_points.json';OUT=BASE/'low_order_segment_masking_points_readback.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def key(r):return (r['network'],tuple(r['segment']),r['location'],r['order'],r['kind'])
def main():
    if OUT.exists():raise FileExistsError(OUT)
    d=json.loads(DATA.read_text())
    if sha(SRC)!=d['source_records_sha256'] or sha(ROOT/'scripts/74_polynomial_forcing_prototype.py')!=d['poly_sha256'] or sha(ROOT/'scripts/259_resolve_segment_masking_20261106.py')!=d['script_sha256']:raise ValueError('hash')
    expected={key(r):r for r in json.loads(SRC.read_text()) if r['kind']!='normal' and F(r['source_derivative'])!=0 and r['contained']};seen=set();counts=Counter()
    for r in d['rows']:
        k=key(r)
        if k not in expected or k in seen:raise ValueError('case')
        seen.add(k);orig=expected[k]
        if (r['fixture_sha256'],r['x'],r['segment_lower'],r['segment_upper'],r['reference110'])!=(orig['fixture_sha256'],orig['x'],orig['lower'],orig['upper'],orig['reference110']):raise ValueError('mismatch')
        lo,hi,ref=map(F,(r['point_lower'],r['point_upper'],r['reference110']));detected=not lo<=ref<=hi
        if lo>hi or not detected or r['detected']!=detected:raise ValueError('detection')
        counts[r['kind']]+=1
    if seen!=set(expected) or len(seen)!=107 or d['records']!=len(seen) or d['detected']!=len(seen) or counts!={'omit_source':59,'reverse_source_sign':48}:raise ValueError('incomplete')
    out={'records':len(seen),'detected':len(seen),'omit_source':counts['omit_source'],'reverse_source_sign':counts['reverse_source_sign'],'records_sha256':sha(DATA),'scope':'readback of selected synthetic cases; numerical point references'}
    with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
    print(json.dumps(out))
if __name__=='__main__':main()
