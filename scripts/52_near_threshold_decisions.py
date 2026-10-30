"""Stress public-CLI decisions at thresholds 1% above/below independent references.

This experiment follows fixture generation and prior CLI observations; it is
not pre-registered. Independent references are used after each CLI completes.
"""
from datetime import datetime, timezone
from fractions import Fraction as F
from decimal import Decimal, localcontext
from pathlib import Path
import hashlib, json, os, subprocess, sys

ROOT=Path(__file__).resolve().parents[1]
FIXTURE=ROOT/'configs/frozen_networks_holdout_2026-10-30.json'

def frac(e):return F(int(e['numerator']),int(e['denominator']))

def main():
    fixture=json.loads(FIXTURE.read_text());gen_dir=ROOT/fixture['config']['generator_run']
    gen=json.loads((gen_dir/'config.json').read_text());refs=json.loads((gen_dir/'network_references.json').read_text())
    if gen['status']!='complete' or hashlib.sha256(FIXTURE.read_bytes()).hexdigest()!=gen['fixture_sha256']:raise ValueError('Invalid fixture generation')
    out=ROOT/'results/generated/52_near_threshold_decisions'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');out.mkdir(parents=True,exist_ok=False)
    template=json.loads((ROOT/'configs/examples/deep_sine.json').read_text());cases=out/'cases';cases.mkdir();records=[]
    config={'status':'running','fixture_sha256':gen['fixture_sha256'],'case_labels':[x['label'] for x in fixture['networks']],
            'degrees':[2,3],'tolerances':'For each case: 0.99 × reference lower decimal, 1.01 × reference upper decimal',
            'scope':'Public CLI decision; thresholds explicitly selected from independent reference after fixture generation. Stress test, not blind holdout. Degree 2/3, cap 4096; no timing claims.',
            'driver_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (out/'executed_source.py').write_bytes(Path(__file__).read_bytes())
    def save():
        (out/'config.json').write_text(json.dumps(config,indent=2)+'\n');(out/'records.json').write_text(json.dumps(records,indent=2)+'\n')
    save()
    try:
        for row in fixture['networks']:
            label=row['label'];case=json.loads(json.dumps(template));case['name']=label;case['problem']['frequency']=row['frequency']
            case['network'].update(fixture=os.path.relpath(FIXTURE,cases).replace('\\','/'),label=label,sha256=gen['fixture_sha256'])
            path=cases/(label+'.json');path.write_text(json.dumps(case,indent=2)+'\n')
            reference=refs[label]
            with localcontext() as ctx:
                ctx.prec=80
                tests=[(str(Decimal(reference['lower']['decimal_outward'])*Decimal('0.99')),4096),
                       (str(Decimal(reference['upper']['decimal_outward'])*Decimal('1.01')),4096)]
            for degree in (2,3):
                for tolerance,cap in tests:
                    output=out/'runs'/(label+'_d'+str(degree)+'_'+tolerance.replace('.','_')+('_cap'+str(cap) if cap else ''))
                    cmd=[sys.executable,str(ROOT/'scripts/certify.py'),'--case',str(path),'--goal','decision','--degree',str(degree),'--error-tolerance',tolerance,'--out',str(output)]
                    cmd += ['--max-intervals',str(cap)]
                    proc=subprocess.run(cmd,cwd=ROOT,text=True,encoding='utf-8',capture_output=True)
                    if proc.returncode not in (0,2):raise RuntimeError(f'CLI failed: {label} degree{degree} {tolerance}: {proc.stderr}')
                    result=json.loads((output/'result.json').read_text());reference=refs[label];lo,hi=frac(reference['lower']),frac(reference['upper']);tau=F(tolerance)
                    decision=result['accuracy_decision'];enclosure=result['enclosure']
                    if enclosure is not None:assert frac(enclosure['lower'])<=lo<=hi<=frac(enclosure['upper'])
                    expected='within_tolerance' if hi<=tau else 'exceeds_tolerance' if lo>tau else 'reference_crosses_threshold'
                    if decision!='undetermined':assert decision==expected,(label,degree,tolerance,decision,expected)
                    assert bool(result['goal_met']) == (decision in ('within_tolerance','exceeds_tolerance'))
                    assert proc.returncode == (0 if result['goal_met'] else 2)
                    records.append(dict(label=label,degree=degree,frequency=row['frequency'],tolerance=tolerance,cap=cap,decision=decision,
                                        reference_comparison=expected,intervals=enclosure['intervals'] if enclosure else None,goal_met=result['goal_met'],
                                        exit_code=proc.returncode,output_dir=str(output.relative_to(ROOT))))
                    save();print('DECISION',label,'degree',degree,tolerance,cap,decision,'segments',records[-1]['intervals'],flush=True)
        config.update(status='complete',requests=len(records),goal_met=sum(r['goal_met'] for r in records));save()
    except Exception as exc:
        config.update(status='failed',error=repr(exc));save();raise
    print('COMPLETED',out,'requests',len(records),flush=True)

if __name__=='__main__':main()
