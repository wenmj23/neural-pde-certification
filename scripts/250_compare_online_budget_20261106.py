"""Post-hoc compare fallback at cap n with unsplit baseline at cap 2n.

Only paired cases where unsplit cap n failed are included. Underlying
derivative calls are counts, not equivalent wall-time or a cost proof.
"""
import hashlib,json
from itertools import product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results/2026-11-06';OUT=B/'fallback_vs_expand_readback.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    if OUT.exists():raise FileExistsError(OUT)
    result={}
    for name,folder,reader in [('polynomial','online_split_fallback','online_split_fallback_readback.json'),('sine','online_sine_fallback','online_sine_fallback_readback.json')]:
        d=B/folder;s=json.loads((d/'summary.json').read_text());q=json.loads((B/reader).read_text());path=d/'records.json'
        if sha(path)!=s['records_sha256'] or sha(path)!=q['records_sha256']:raise ValueError('records hash')
        rows=json.loads(path.read_text());idx={(r['network'],r['degree'],r['cap']):r for r in rows};names={r['network'] for r in rows};groups={}
        for degree,cap in product((2,3),(32,64)):
            pairs=[]
            for n in names:
                low,high=idx[n,degree,cap],idx[n,degree,cap*2]
                if not low['fallback_called']:continue
                low_base=[a for a in low['attempts'] if a['pieces']==1];low_extra=[a for a in low['attempts'] if a['pieces']==2];high_base=[a for a in high['attempts'] if a['pieces']==1]
                if not low_extra or low_base[-1]['goal_met']:raise ValueError('unexpected fallback')
                pairs.append({'network':n,'fallback_success':low['goal_met'],'expanded_base_success':high_base[-1]['goal_met'],
                              'fallback_calls':low_base[-1]['derivative_calls_cumulative']+low_extra[-1]['derivative_calls_cumulative'],
                              'expanded_base_calls':high_base[-1]['derivative_calls_cumulative']})
            groups[f'{degree}/{cap}_to_{2*cap}']={'cases':len(pairs),'both_success':sum(p['fallback_success'] and p['expanded_base_success'] for p in pairs),
                'only_fallback':sum(p['fallback_success'] and not p['expanded_base_success'] for p in pairs),'only_expanded':sum(not p['fallback_success'] and p['expanded_base_success'] for p in pairs),
                'neither':sum(not p['fallback_success'] and not p['expanded_base_success'] for p in pairs),'fallback_more_calls':sum(p['fallback_calls']>p['expanded_base_calls'] for p in pairs),'pairs':pairs}
        result[name]={'groups':groups,'source_records_sha256':sha(path)}
    out={'datasets':result,'scope':'post-hoc comparison after targets and network outcomes known; derivative call count only, not matched time or general blind performance'}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':main()
