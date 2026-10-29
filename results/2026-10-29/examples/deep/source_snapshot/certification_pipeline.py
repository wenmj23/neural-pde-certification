"""Reproducible CLI workflow around the existing scoped Poisson certifiers.

No training, exact solution, reference file or observed error is used to select
a mesh or stop. All decisions compare exact rational interval endpoints.
"""
from copy import deepcopy
from datetime import datetime, timezone
from decimal import Decimal, localcontext, ROUND_CEILING
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import re
import time

import flint
import mpmath
from arb_certification import export
from arb_cubic_reconstruction import CubicVerifier, SixthOrderDeep
from arb_shaped_remainder import ShapedVerifier

ROOT = Path(__file__).resolve().parents[1]
SOURCE_NAMES = [
    'certify.py', 'certification_pipeline.py', 'arb_cubic_reconstruction.py',
    'arb_deep_network.py', 'arb_shaped_remainder.py', 'arb_reconstruction.py',
    'arb_certification.py', 'arb_problem_families.py', 'interval_certification.py',
]
DEFAULTS = dict(degree=2, mesh='uniform', goal='enclosure', gap_ratio='1.02',
                error_tolerance=None, initial_intervals=16, max_intervals=4096,
                max_attempts=10, time_limit_seconds='30', precision_bits=103)
SCOPE = [
    'Unit interval, -u_dd=(k*pi)^2*sin(k*pi*x), zero Dirichlet data.',
    'Exact mathematical dense tanh network with frozen binary64 weights and output x*(1-x)*N(x).',
    'Quadratic requires residual C3; cubic requires residual C4. These tanh/sine cases are smooth.',
    'Bounds are conditional on the stated mathematics and interval backend inclusion correctness; no independent/formal software audit.',
]


def digest(data):
    return hashlib.sha256(data).hexdigest()


def strict_json(data):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError('Duplicate JSON key: ' + key)
            result[key] = value
        return result
    def invalid(value):
        raise ValueError('Nonfinite JSON number: ' + value)
    return json.loads(data, object_pairs_hook=pairs, parse_constant=invalid)


def fields(value, allowed, required, name):
    if not isinstance(value, dict):
        raise ValueError(name + ' must be an object')
    extra = set(value) - set(allowed)
    missing = set(required) - set(value)
    if extra or missing:
        raise ValueError(f'{name}: unknown fields {sorted(extra)}, missing fields {sorted(missing)}')


def rational(value, name, minimum=Fraction(0), strict=True):
    # Floats would silently change the user's requested threshold.
    if not isinstance(value, str) or not value or len(value) > 120:
        raise ValueError(name + ' must be a finite decimal or fraction string')
    try:
        result = Fraction(value)
    except (ValueError, ZeroDivisionError, OverflowError) as error:
        raise ValueError('Invalid ' + name) from error
    if result < minimum or (strict and result == minimum):
        raise ValueError(name + ' is below the allowed minimum')
    return result


def integer(value, name, lower, upper):
    if type(value) is not int or not lower <= value <= upper:
        raise ValueError(f'{name} must be an integer in [{lower}, {upper}]')
    return value


def validate_settings(settings):
    fields(settings, DEFAULTS, (), 'verification')
    result = {**DEFAULTS, **settings}
    integer(result['degree'], 'degree', 2, 3)
    if result['mesh'] not in ('uniform', 'adaptive'):
        raise ValueError('mesh must be uniform or adaptive')
    if result['goal'] not in ('enclosure', 'decision'):
        raise ValueError('goal must be enclosure or decision')
    rational(result['gap_ratio'], 'gap_ratio', Fraction(1))
    if result['error_tolerance'] is not None:
        rational(result['error_tolerance'], 'error_tolerance', strict=False)
    if result['goal'] == 'decision' and result['error_tolerance'] is None:
        raise ValueError('decision goal requires error_tolerance')
    integer(result['initial_intervals'], 'initial_intervals', 1, 16384)
    integer(result['max_intervals'], 'max_intervals', result['initial_intervals'], 16384)
    integer(result['max_attempts'], 'max_attempts', 1, 32)
    integer(result['precision_bits'], 'precision_bits', 80, 700)
    if rational(result['time_limit_seconds'], 'time_limit_seconds') > 3600:
        raise ValueError('time_limit_seconds must not exceed 3600')
    return result


def load_case(path, overrides=None):
    path = Path(path).resolve()
    original = path.read_bytes()
    case = strict_json(original)
    fields(case, ('schema_version', 'name', 'problem', 'network', 'verification'),
           ('schema_version', 'name', 'problem', 'network'), 'case')
    if type(case['schema_version']) is not int or case['schema_version'] != 1:
        raise ValueError('Unsupported case schema_version')
    if (not isinstance(case['name'], str)
            or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,79}', case['name'])):
        raise ValueError('Invalid case name')
    p = case['problem']
    fields(p, ('operator', 'source', 'frequency', 'domain', 'boundary'),
           ('operator', 'source', 'frequency', 'domain', 'boundary'), 'problem')
    if p['operator'] != '-u_dd' or p['source'] != '(k*pi)^2*sin(k*pi*x)':
        raise ValueError('Only the sine-forced 1D Poisson family is supported')
    if (p['domain'] != [0, 1] or p['boundary'] != [0, 0]
            or any(type(x) is not int for x in p['domain'] + p['boundary'])):
        raise ValueError('Exactly domain [0,1] and zero Dirichlet boundary [0,0] required')
    integer(p['frequency'], 'frequency', 1, 10000)
    n = case['network']
    fields(n, ('fixture', 'label', 'sha256', 'output_transform'),
           ('fixture', 'label', 'sha256', 'output_transform'), 'network')
    if n['output_transform'] != 'x*(1-x)*raw_network(x)':
        raise ValueError('Unsupported output transform')
    if not isinstance(n['fixture'], str) or not n['fixture'] or not isinstance(n['label'], str):
        raise ValueError('fixture and label must be strings')
    if not isinstance(n['sha256'], str) or not re.fullmatch(r'[a-f0-9]{64}', n['sha256']):
        raise ValueError('Expected a pinned SHA-256 for the frozen fixture')
    fixture = (path.parent / n['fixture']).resolve()
    fixture_bytes = fixture.read_bytes()
    if digest(fixture_bytes) != n['sha256']:
        raise ValueError('Frozen fixture checksum mismatch')
    data = strict_json(fixture_bytes)
    if not isinstance(data, dict) or not isinstance(data.get('networks'), list):
        raise ValueError('Invalid dense-network fixture')
    matches = [row for row in data['networks'] if isinstance(row, dict) and row.get('label') == n['label']]
    if len(matches) != 1:
        raise ValueError('Network label is missing or ambiguous')
    row = matches[0]
    if type(row.get('frequency')) is not int or row['frequency'] != p['frequency']:
        raise ValueError('Case PDE frequency differs from the frozen network problem')
    given = case.get('verification', {})
    fields(given, DEFAULTS, (), 'verification')
    if overrides is not None:
        fields(overrides, DEFAULTS, (), 'overrides')
    settings = validate_settings({**given, **(overrides or {})})
    case = deepcopy(case)
    case['verification'] = settings
    return dict(path=path, original=original, case=case, fixture=fixture,
                fixture_bytes=fixture_bytes, row=row)


def endpoint(record):
    return Fraction(int(record['numerator']), int(record['denominator']))


def classify(lower, upper, settings):
    if lower < 0 or upper < lower:
        raise ArithmeticError('Invalid error-norm enclosure')
    tolerance = settings['error_tolerance']
    if tolerance is None:
        decision = 'not_requested'
    else:
        tolerance = rational(tolerance, 'error_tolerance', strict=False)
        decision = ('within_tolerance' if upper <= tolerance else
                    'exceeds_tolerance' if lower > tolerance else 'undetermined')
    gap_met = ((lower == upper == 0)
               or (lower > 0 and upper <= rational(settings['gap_ratio'], 'gap_ratio', Fraction(1)) * lower))
    goal_met = gap_met if settings['goal'] == 'enclosure' else decision in ('within_tolerance', 'exceeds_tolerance')
    return dict(accuracy_decision=decision, gap_target_met=bool(gap_met),
                goal_met=bool(goal_met), gap_ratio_upper=str(upper/lower) if lower else None)


def serialize(r):
    return dict(lower=export(r['lower'], 'down'), upper=export(r['upper'], 'up'),
                z_lower=export(r['z_lower'], 'down'), z_upper=export(r['z_upper'], 'up'),
                remainder_upper=export(r['remainder_upper'], 'up'), intervals=r['intervals'],
                residual_evaluations=r['residual_evaluations'],
                derivative_interval_evaluations=r['curvature_evaluations'])


def human_summary(result):
    final = result['enclosure']
    decisions = dict(within_tolerance='误差已证明不超过容差', exceeds_tolerance='误差已证明超过容差',
                     undetermined='当前范围尚不能判断是否满足容差', not_requested='未设置绝对误差容差')
    if final is None:
        interval = '未得到有限误差范围'
    else:
        interval = f"最大绝对误差 ∈ [{final['lower']['decimal_outward']}, {final['upper']['decimal_outward']}]"
    settings = result['settings']
    goal = ('收紧误差范围，上界/下界≤' + settings['gap_ratio']
            if settings['goal'] == 'enclosure' else '判断误差是否不超过容差' + settings['error_tolerance'])
    reasons = dict(mesh_limit='已用尽网格预算', time_limit='已用尽时间预算（软限制）',
                   attempt_limit='已用尽尝试次数', goal_met='本次目标已达到')
    return ['本次目标：' + goal, interval, decisions[result['accuracy_decision']],
            ('本次目标已达到' if result['goal_met'] else '本次目标未达到：' + reasons[result['stop_reason']])]


def write_report(output, result):
    settings = result['settings']
    final = result['enclosure']
    lines = ['# 一维 Poisson 最大误差验证结果', '', *human_summary(result), '',
             f"- 案例：{result['case_name']}；网络：{result['network']['label']}",
             f"- 方程：−u''=(kπ)²sin(kπx)，k={result['problem']['frequency']}；区域[0,1]；两端为零",
             '- 近似解：v(x)=x(1−x)N(x)，N为冻结的tanh网络',
             f"- 方法：{settings['degree']}次残差重构＋形状余项包络；网格：{settings['mesh']}",
             f"- 目标：{settings['goal']}；上下界比例要求：{settings['gap_ratio']}",
             f"- 绝对最大误差容差：{settings['error_tolerance']}",
             f"- 终止原因：{result['stop_reason']}；尝试次数：{result['attempts']}", '',
             '上下界比例衡量误差范围有多窄，不表示网络解本身的相对误差。',
             '已证明超过容差是对解的否定判断；预算用尽且范围跨越容差表示尚不能判定。', '',
             '## 计算过程', '',
             '| 尝试 | 区间段数 | 连续误差下界 | 连续误差上界 | 状态 |',
             '|---|---:|---:|---:|---|']
    for row in result['trace']:
        if row['status'] == 'finite':
            lines.append(f"| {row['attempt']} | {row['intervals']} | {row['enclosure']['lower']['decimal_outward']} | {row['enclosure']['upper']['decimal_outward']} | {row['accuracy_decision']} |")
        else:
            lines.append(f"| {row['attempt']} | {row['intervals']} | — | — | 区间运算未完成，未据此作决定 |")
    if final:
        lines += ['', '## 最后一个完整范围的来源', '',
                  f"- 重构函数范数下界：{final['z_lower']['decimal_outward']}",
                  f"- 重构函数范数上界：{final['z_upper']['decimal_outward']}",
                  f"- 重构余项传播上界：{final['remainder_upper']['decimal_outward']}",
                  '', '最终范围为 [max(0, Z_L−D_U), Z_U+D_U]。',
                  f"该范围对应{final['intervals']}段网格；完整坐标保存在mesh.json。"]
    lines += ['', '## 条件与复现', '',
              '这是连续区域上的条件性区间保证；依赖数学假设和区间后端的包含正确性，尚无独立或形式化软件审计。',
              '验证对象是冻结权重定义的精确数学网络。训练、真实解和密集网格观察值不参与判定。',
              '时间预算在完整网格尝试之间检查，最后一次尝试可能超过该软上限。', '',
              'result.json保留精确有理端点、版本和源码摘要；trace.json保留所有尝试。',
              'case.json与fixture.json是本次有效输入；case.original.json保留原始配置。',
              '用统一入口重新运行case.json即可复现；source_snapshot保存本次执行源码。', '']
    (output / 'report.md').write_text('\n'.join(lines), encoding='utf-8')


def run_case(case_path, output=None, overrides=None):
    started = time.perf_counter()
    loaded = load_case(case_path, overrides)
    case = loaded['case']; settings = case['verification']; network = case['network']
    # Reuse the reviewed research kernels without changing historical sources.
    problem = SixthOrderDeep(loaded['fixture'], network['label'], network['sha256'], settings['precision_bits'])
    verifier = (ShapedVerifier(problem, 'moment') if settings['degree'] == 2 else CubicVerifier(problem, 'moment'))
    if output is None:
        output = ROOT / 'results/generated/certify' / case['name'] / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output / 'case.original.json').write_bytes(loaded['original'])
    replay = deepcopy(case); replay['network']['fixture'] = 'fixture.json'
    replay_bytes = (json.dumps(replay, indent=2, ensure_ascii=False) + '\n').encode()
    (output / 'case.json').write_bytes(replay_bytes)
    (output / 'fixture.json').write_bytes(loaded['fixture_bytes'])
    snapshots = output / 'source_snapshot'; snapshots.mkdir()
    sources = {}
    for name in SOURCE_NAMES:
        source = Path(__file__).resolve().parent / name
        content = source.read_bytes(); sources[name] = digest(content)
        (snapshots / name).write_bytes(content)
    n = settings['initial_intervals']; verifier.uniform(n)
    trace = []; final = None; final_nodes = None; stop = 'attempt_limit'
    classification = dict(accuracy_decision='undetermined' if settings['error_tolerance'] is not None else 'not_requested',
                          gap_target_met=False, goal_met=False, gap_ratio_upper=None)
    deadline = float(rational(settings['time_limit_seconds'], 'time_limit_seconds'))
    for attempt in range(1, settings['max_attempts'] + 1):
        if time.perf_counter() - started >= deadline:
            stop = 'time_limit'; break
        attempt_started = time.perf_counter()
        row = dict(attempt=attempt, intervals=n)
        try:
            if settings['mesh'] == 'adaptive':
                verifier.refine_to(n)
            else:
                verifier.uniform(n)
            current = serialize(verifier.certify())
            classification = classify(endpoint(current['lower']), endpoint(current['upper']), settings)
            final = current; final_nodes = list(verifier.nodes)
            row.update(status='finite', enclosure=current, **classification)
        except ArithmeticError as error:
            row.update(status='arithmetic_unresolved', reason=str(error))
        row['elapsed_seconds'] = time.perf_counter() - attempt_started
        trace.append(row)
        # Complete snapshots plus readback avoid silently losing the last checks.
        (output / 'trace.json').write_text(json.dumps(trace, indent=2) + '\n')
        if row['status'] == 'finite' and classification['goal_met']:
            stop = 'goal_met'; break
        if time.perf_counter() - started >= deadline:
            stop = 'time_limit'; break
        if n == settings['max_intervals']:
            stop = 'mesh_limit'; break
        n = min(2*n, settings['max_intervals'])
    result = dict(schema_version=1, case_name=case['name'], error_norm='Linf_absolute',
                  problem=case['problem'], settings=settings,
                  network=dict(label=network['label'], hidden_widths=[len(x['bias']) for x in loaded['row']['layers'][:-1]],
                               fixture_sha256=network['sha256'], output_transform=network['output_transform']),
                  stop_reason=stop, enclosure=final, attempts=len(trace), **classification,
                  work=dict(residual_evaluations=verifier.residual_evaluations,
                            derivative_interval_evaluations=verifier.curvature_evaluations),
                  wall_seconds=time.perf_counter()-started,
                  versions=dict(python=platform.python_version(), python_flint=flint.__version__,
                                flint=flint.__FLINT_VERSION__, mpmath=mpmath.__version__),
                  input_sha256=dict(original_case=digest(loaded['original']), replay_case=digest(replay_bytes)),
                  source_sha256=sources, assumptions=SCOPE, trace=trace)
    (output / 'mesh.json').write_text(json.dumps(dict(nodes=[str(x) for x in final_nodes] if final_nodes else None), indent=2)+'\n')
    (output / 'trace.json').write_text(json.dumps(trace, indent=2)+'\n')
    (output / 'result.json').write_text(json.dumps(result, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    write_report(output, result)
    saved = strict_json((output / 'result.json').read_bytes())
    saved_trace = strict_json((output / 'trace.json').read_bytes())
    if saved != result or saved_trace != trace:
        raise OSError('Output readback mismatch')
    return result, output
