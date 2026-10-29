"""Outcome-oriented checks for the public certification workflow.

Exercise exact threshold decisions, inconclusive budgets, input/PDE integrity,
reference independence, CLI exit codes and replay from a portable snapshot.
"""
from copy import deepcopy
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import certification_pipeline as workflow

ROOT = Path(__file__).resolve().parents[1]
SHALLOW = ROOT / 'configs/examples/shallow_sine.json'
DEEP = ROOT / 'configs/examples/deep_sine.json'


def main():
    out = ROOT / 'results/generated/32_pipeline_checks' / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True, exist_ok=False)
    checks = []
    sources = [Path(__file__)] + [ROOT/'scripts'/name for name in workflow.SOURCE_NAMES]
    source_hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    snapshot = out/'test_source_snapshot'; snapshot.mkdir()
    for p in sources:
        (snapshot/p.name).write_bytes(p.read_bytes())
    def save(name, value):
        (out/name).write_text(json.dumps(value, indent=2, ensure_ascii=False)+'\n')
    save('config.json', dict(status='running', sources=source_hashes))
    def check(name, passed):
        checks.append(dict(name=name, passed=bool(passed)))
        save('checks.json', checks)
        if not passed:
            raise AssertionError(name)
    # Exact comparisons include the equality case at the user's tolerance.
    settings = workflow.validate_settings(dict(goal='decision', error_tolerance='1/1000'))
    for name, lo, hi, decision in [
        ('upper_equals_tolerance', F(0), F(1,1000), 'within_tolerance'),
        ('lower_equals_tolerance', F(1,1000), F(2,1000), 'undetermined'),
        ('strictly_above', F(2,1000), F(3,1000), 'exceeds_tolerance'),
        ('straddles', F(1,2000), F(3,2000), 'undetermined'),
    ]:
        r = workflow.classify(lo, hi, settings)
        check(name, r['accuracy_decision'] == decision)
        check(name+'_goal', r['goal_met'] == (decision in ('within_tolerance','exceeds_tolerance')))
    settings = workflow.validate_settings({})
    check('zero_exact_gap', workflow.classify(F(0), F(0), settings)['gap_target_met'])
    check('zero_lower_not_tight', not workflow.classify(F(0), F(1), settings)['gap_target_met'])
    check('gap_equality', workflow.classify(F(1), F(51,50), settings)['gap_target_met'])
    invalid = out/'invalid'; invalid.mkdir()
    baseline = json.loads(SHALLOW.read_text())
    baseline['network']['fixture'] = str((SHALLOW.parent/baseline['network']['fixture']).resolve())
    mutations = [
        ('wrong_pde_frequency', lambda c: c['problem'].__setitem__('frequency', 5)),
        ('wrong_source', lambda c: c['problem'].__setitem__('source', '2')),
        ('wrong_boundary', lambda c: c['problem'].__setitem__('boundary', [1,0])),
        ('wrong_domain', lambda c: c['problem'].__setitem__('domain', [0,2])),
        ('wrong_checksum', lambda c: c['network'].__setitem__('sha256', '0'*64)),
        ('unknown_label', lambda c: c['network'].__setitem__('label', 'missing')),
        ('wrong_transform', lambda c: c['network'].__setitem__('output_transform', 'raw_network(x)')),
        ('boolean_schema', lambda c: c.__setitem__('schema_version', True)),
        ('verification_not_object', lambda c: c.__setitem__('verification', [])),
        ('unknown_setting', lambda c: c['verification'].__setitem__('max_intervlas', 16)),
        ('float_threshold', lambda c: c['verification'].__setitem__('gap_ratio', 1.02)),
        ('nonfinite_threshold', lambda c: c['verification'].__setitem__('error_tolerance', 'nan')),
        ('negative_threshold', lambda c: c['verification'].__setitem__('error_tolerance', '-1')),
        ('zero_denominator', lambda c: c['verification'].__setitem__('gap_ratio', '1/0')),
        ('impossible_gap', lambda c: c['verification'].__setitem__('gap_ratio', '1')),
        ('invalid_degree', lambda c: c['verification'].__setitem__('degree', 4)),
        ('invalid_precision', lambda c: c['verification'].__setitem__('precision_bits', 20)),
        ('boolean_attempts', lambda c: c['verification'].__setitem__('max_attempts', True)),
        ('small_mesh_limit', lambda c: c['verification'].__setitem__('max_intervals', 8)),
        ('missing_decision_tolerance', lambda c: c['verification'].update(goal='decision', error_tolerance=None)),
    ]
    for name, mutate in mutations:
        case = deepcopy(baseline); mutate(case); path=invalid/(name+'.json'); path.write_text(json.dumps(case))
        rejected = False
        try:
            workflow.load_case(path)
        except (ValueError, TypeError):
            rejected = True
        check('reject_'+name, rejected)
    for name, data in [('duplicate_key', '{"a":1,"a":2}'), ('json_nan', '{"a":NaN}')]:
        rejected=False
        try: workflow.strict_json(data)
        except ValueError: rejected=True
        check(name, rejected)
    reference_path = ROOT/'results/generated/24_depth_reference/20261029T044806546934Z/references.json'
    references = json.loads(reference_path.read_text())
    sample_results = {}
    for name, path in [('shallow',SHALLOW), ('deep',DEEP)]:
        r, folder=workflow.run_case(path, out/name); sample_results[name]=(r,folder)
        ref=references[r['network']['label']]; final=r['enclosure']
        check(name+'_goal_met', r['goal_met'] and r['stop_reason']=='goal_met')
        check(name+'_continuous_reference_contained', workflow.endpoint(final['lower'])<=workflow.endpoint(ref['lower']) and workflow.endpoint(final['upper'])>=workflow.endpoint(ref['upper']))
        check(name+'_accuracy_decision', r['accuracy_decision']==('within_tolerance' if name=='shallow' else 'exceeds_tolerance'))
        check(name+'_snapshots', all(hashlib.sha256((folder/'source_snapshot'/n).read_bytes()).hexdigest()==h for n,h in r['source_sha256'].items()))
        check(name+'_last_mesh_matches', len(json.loads((folder/'mesh.json').read_text())['nodes'])==final['intervals']+1)
        print('WORKFLOW',name,r['accuracy_decision'],final['intervals'],flush=True)
    for degree in (2,3):
        r,_=workflow.run_case(SHALLOW,out/f'adaptive_degree{degree}',dict(degree=degree,mesh='adaptive'))
        check('adaptive_goal_'+str(degree),r['goal_met'] and r['accuracy_decision']=='within_tolerance')
    r,_=workflow.run_case(DEEP,out/'quadratic_deep',dict(degree=2))
    check('quadratic_deep',r['goal_met'] and r['accuracy_decision']=='exceeds_tolerance')
    r,_=workflow.run_case(SHALLOW,out/'mesh_limit',dict(max_intervals=16,error_tolerance='0.0001'))
    check('mesh_limit_is_unknown',not r['goal_met'] and r['stop_reason']=='mesh_limit' and r['accuracy_decision']=='undetermined')
    r,_=workflow.run_case(SHALLOW,out/'attempt_limit',dict(max_attempts=1,error_tolerance='0.0001'))
    check('attempt_limit_is_unknown',not r['goal_met'] and r['stop_reason']=='attempt_limit' and r['accuracy_decision']=='undetermined')
    r,_=workflow.run_case(SHALLOW,out/'time_limit',dict(time_limit_seconds='1e-30'))
    check('no_attempt_is_not_zero_error',not r['goal_met'] and r['stop_reason']=='time_limit' and r['enclosure'] is None and r['accuracy_decision']=='undetermined')
    r,_=workflow.run_case(SHALLOW,out/'decision_shallow',dict(goal='decision'))
    check('decision_can_precede_tightness',r['goal_met'] and not r['gap_target_met'] and r['enclosure']['intervals']<sample_results['shallow'][0]['enclosure']['intervals'])
    r,_=workflow.run_case(DEEP,out/'decision_deep',dict(goal='decision'))
    check('negative_decision_is_goal_met',r['goal_met'] and r['accuracy_decision']=='exceeds_tolerance')
    # Minimal fixture has weights and PDE label only, no observed/true errors.
    loaded=workflow.load_case(SHALLOW); original_row=loaded['row']
    minimal=dict(networks=[{k:original_row[k] for k in ('label','frequency','layers')}])
    minimal_path=out/'minimal_fixture.json';minimal_path.write_text(json.dumps(minimal))
    case=deepcopy(baseline);case['network'].update(fixture=str(minimal_path),sha256=hashlib.sha256(minimal_path.read_bytes()).hexdigest())
    minimal_case=out/'minimal_case.json';minimal_case.write_text(json.dumps(case))
    r,_=workflow.run_case(minimal_case,out/'no_reference_fields')
    check('no_reference_fields_needed',r['enclosure']['lower']==sample_results['shallow'][0]['enclosure']['lower'] and r['enclosure']['upper']==sample_results['shallow'][0]['enclosure']['upper'])
    # Existing output directories are never overwritten.
    preserved=(out/'shallow/result.json').read_bytes();rejected=False
    try:workflow.run_case(SHALLOW,out/'shallow')
    except FileExistsError:rejected=True
    check('no_overwrite',rejected and (out/'shallow/result.json').read_bytes()==preserved)
    # Portable replay uses only the copied inputs and copied source tree.
    portable=out/'portable_bundle';shutil.copytree(sample_results['shallow'][1],portable)
    empty=out/'empty_cwd';empty.mkdir()
    command=[sys.executable,str(portable/'source_snapshot/certify.py'),'--case',str(portable/'case.json'),'--out',str(out/'portable_replay')]
    process=subprocess.run(command,cwd=empty,text=True,capture_output=True,env=os.environ.copy())
    (out/'portable_stdout.txt').write_text(process.stdout);(out/'portable_stderr.txt').write_text(process.stderr)
    check('portable_cli_exit',process.returncode==0)
    replay=json.loads((out/'portable_replay/result.json').read_text())
    check('portable_exact_endpoints',replay['enclosure']['lower']==sample_results['shallow'][0]['enclosure']['lower'] and replay['enclosure']['upper']==sample_results['shallow'][0]['enclosure']['upper'])
    for name, args, expected in [
        ('decision_accept',['--case',str(SHALLOW),'--goal','decision'],0),
        ('decision_reject',['--case',str(DEEP),'--goal','decision'],0),
        ('unknown',['--case',str(SHALLOW),'--max-intervals','16','--error-tolerance','0.0001'],2),
        ('invalid',['--case',str(invalid/'wrong_pde_frequency.json')],1),
    ]:
        command=[sys.executable,str(ROOT/'scripts/certify.py'),*args,'--out',str(out/('cli_'+name))]
        process=subprocess.run(command,cwd=empty,text=True,capture_output=True,env=os.environ.copy())
        (out/('cli_'+name+'.txt')).write_text(process.stdout+'\n'+process.stderr)
        check('cli_exit_'+name,process.returncode==expected)
    check('no_training_or_plot_imports',not any(n in sys.modules for n in ('torch','numpy','matplotlib')))
    save('checks.json',checks)
    saved=json.loads((out/'checks.json').read_text());check_count=len(checks)
    if len(saved)!=check_count or not all(x['passed'] for x in saved):raise AssertionError('Check file incomplete')
    config=dict(status='complete',checks_passed=check_count,sources=source_hashes,
                examples={name:str(folder.relative_to(ROOT)) for name,(_,folder) in sample_results.items()},
                continuous_reference=str(reference_path.relative_to(ROOT)),
                scope='Workflow checks around existing conditional interval kernels; no new PDE/architecture theorem or independent software audit')
    save('config.json',config)
    print('COMPLETED',out,'checks',check_count,flush=True)


if __name__=='__main__':main()
