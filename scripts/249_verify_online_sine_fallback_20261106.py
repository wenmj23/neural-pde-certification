"""Read back online fallback traces and compare prior fixed strategies."""
from fractions import Fraction as F
import hashlib,json,math
from itertools import product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results/2026-11-06';D=B/'online_sine_fallback';OUT=B/'online_sine_fallback_readback.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(x,s):
    if not x:raise ValueError(s)
def val(x):return F(int(x['numerator']),int(x['denominator']))
def main():
    if OUT.exists():raise FileExistsError(OUT)
    c=json.loads((D/'locked_config.json').read_text());s=json.loads((D/'summary.json').read_text());rows=json.loads((D/'records.json').read_text());prior=ROOT/c['prior'];ic=ROOT/c['input'];inputs=json.loads(ic.read_text())
    need(sha(D/'records.json')==s['records_sha256'] and sha(D/'locked_config.json')==s['config_sha256']==sha(ROOT/'configs/online_sine_fallback_2026-11-06.json') and sha(ROOT/'scripts/248_online_sine_fallback_20261106.py')==s['script_sha256'],'data/code')
    need(sha(prior/'records.json')==s['prior_records_sha256'] and sha(prior/'locked_config.json')==sha(ic) and sha(ROOT/'scripts/225_derivative_split_20261105.py')==s['split_source_sha256'],'source')
    f=ROOT/inputs['fixture'];rp=ROOT/inputs['reference'];need(sha(f)==inputs['fixture_sha256'] and sha(rp)==inputs['reference_sha256'],'fixture/reference')
    refs={q['network']:q for q in json.loads(rp.read_text())};nets=json.loads(f.read_text())['networks'];need(len(nets)==12 and set(n['label'] for n in nets)==set(refs),'manifest')
    cases={n['label']:(F(refs[n['label']]['lower']),F(refs[n['label']]['upper'])) for n in nets}
    need(all(refs[n]['converged'] for n in cases),'reference convergence')
    keys=[(r['network'],r['degree'],r['cap']) for r in rows];expected=set(product(cases,c['degrees'],c['caps']));need(len(cases)==12 and len(rows)==len(set(keys))==len(expected)==72 and set(keys)==expected,'keys');old={(r['network'],r['degree'],r['pieces']):r for r in json.loads((prior/'records.json').read_text())};groups={};failures=[];attempt_count=0
    for row in rows:
        need(row['status']=='complete' and row['attempts'],'state');name,d,cap=row['network'],row['degree'],row['cap'];rl,ru=cases[name];allowed=[n for n in c['mesh_sequence'] if n<=cap];base=[a for a in row['attempts'] if a['pieces']==1];extra=[a for a in row['attempts'] if a['pieces']==2]
        need([a['mesh'] for a in base]==allowed[:len(base)] and 1<=len(base)<=len(allowed),'baseline path');need(not any(a['goal_met'] for a in base[:-1]),'baseline stop');need(bool(extra)==(not base[-1]['goal_met'])==row['fallback_called'],'fallback trigger')
        if extra:need(len(base)==len(allowed) and [a['mesh'] for a in extra]==allowed[:len(extra)] and 1<=len(extra)<=len(allowed) and not any(a['goal_met'] for a in extra[:-1]),'fallback path')
        need(row['attempts']==base+extra and row['goal_met']==row['attempts'][-1]['goal_met'],'outcome');need(row['goal_met'] or len(extra)==len(allowed),'premature stop')
        elapsed=0;evals={1:0,2:0}
        for a in row['attempts']:
            attempt_count+=1;need(math.isfinite(a['seconds']) and a['seconds']>0,'seconds');elapsed+=a['seconds'];need(abs(elapsed-a['cumulative_seconds'])<1e-9,'time');evals[a['pieces']]+=a['mesh'];need(a['derivative_calls_cumulative']==evals[a['pieces']]*a['pieces'],'derivative calls')
            lo,hi=F(a['lower']),F(a['upper']);need(0<=lo<=rl<=ru<=hi and a['goal_met']==(lo>0 and hi<=F(c['gap_ratio'])*lo),'bounds/goal');previous=old[(name,d,a['pieces'])]['attempts'];matches=[x for x in previous if x['mesh']==a['mesh']];need(len(matches)==1 and (lo,hi)==(F(matches[0]['lower']),F(matches[0]['upper'])),'prior endpoint')
        oldbase=old[(name,d,1)]['attempts'];oldsplit=old[(name,d,2)]['attempts'];expected_base=any(x['gap_met'] and x['mesh']<=cap for x in oldbase);expected_split=any(x['gap_met'] and x['mesh']<=cap for x in oldsplit)
        need(base[-1]['goal_met']==expected_base and row['goal_met']==(expected_base or expected_split),'simulation vs actual')
        key=f'{d}/{cap}';g=groups.setdefault(key,{'cases':0,'baseline_success':0,'fallback_called':0,'fallback_success':0,'final_success':0,'final_unknown':0});g['cases']+=1;g['baseline_success']+=expected_base;g['fallback_called']+=bool(extra);g['fallback_success']+=bool(extra) and row['goal_met'];g['final_success']+=row['goal_met'];g['final_unknown']+=not row['goal_met']
        if not row['goal_met']:failures.append({'network':name,'degree':d,'cap':cap})
    need(attempt_count==s['attempts'] and sum(r['fallback_called'] for r in rows)==s['fallback_calls'] and sum(r['goal_met'] for r in rows)==s['goal_met'],'summary')
    o={'records':len(rows),'attempts':attempt_count,'matched_prior_paths':len(rows),'groups':groups,'failures':failures,'records_sha256':sha(D/'records.json'),'scope':'post-exit audit of stored rational endpoints, provenance and online branching; not mathematical proof'};OUT.write_text(json.dumps(o,indent=2)+'\n');print(json.dumps(o),flush=True)
if __name__=='__main__':main()
