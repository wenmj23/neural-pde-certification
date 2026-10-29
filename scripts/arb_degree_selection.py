"""Research-only degree-2/3 scheduling over unchanged validated kernels.

Forecasts use measured work and certificate Z/D components, never true errors.
They affect efficiency only: every successful result passes an exact enclosure
test. The public certify.py interface is intentionally unchanged.
"""
from fractions import Fraction as F
import math
import time

from arb_certification import exact_fraction, export
from arb_cubic_reconstruction import CubicVerifier
from arb_shaped_remainder import ShapedVerifier

POLICIES = ('fixed2', 'fixed3', 'probe_cubic', 'forecast_once', 'forecast_race', 'forecast_guarded',
            'microprobe_cubic', 'fixed2_from1', 'fixed3_from1')


def pack(r):
    result = {name: export(r[name], 'down' if name in ('lower', 'z_lower') else 'up')
              for name in ('lower', 'upper', 'z_lower', 'z_upper', 'remainder_upper')}
    result.update(intervals=r['intervals'], degree=r['degree'])
    return result


def fraction(endpoint):
    return F(int(endpoint['numerator']), int(endpoint['denominator']))


def gap_met(r, target):
    lo, hi = exact_fraction(r['lower']), exact_fraction(r['upper'])
    if lo < 0 or hi < lo:
        raise ArithmeticError('Invalid enclosure')
    return lo == hi == 0 or (lo > 0 and hi <= target*lo)


def log_fraction(value):
    if value < 0:
        raise ValueError('Nonnegative forecast quantity required')
    return math.log2(value.numerator)-math.log2(value.denominator) if value else -math.inf


def components(r):
    zl, zu = exact_fraction(r['z_lower']), exact_fraction(r['z_upper'])
    du = exact_fraction(r['remainder_upper'])
    if zl < 0 or zu < zl or du < 0:
        raise ArithmeticError('Invalid forecast components')
    return (zl+zu)/2, (zu-zl)/2, du


def forecast(state, target, cap, observed=False):
    """Predict *additional* cost; this is not a mathematical error bound.

    For a symmetric bracket, gap<=target iff halfwidth/center <=
    (target-1)/(target+1). Predict Z halfwidth as h^2 and remainder as
    h^(degree+1). Optional measured contraction is clamped to [1,12].
    A zero lower bound does not make the predictor falsely report success.
    """
    n = state['n']
    if n >= cap:
        return dict(seconds=None, steps=None, reason='mesh_limit')
    last = state['last']
    elapsed = max(state['last_seconds'], 1e-6)
    if last is None:
        return dict(seconds=2*elapsed, steps=1, reason='arithmetic_retry')
    center, zwidth, remainder = components(last)
    if not center:
        return dict(seconds=None, steps=None, reason='zero_reconstructed_scale')
    zorder, dorder = 2., float(state['degree']+1)
    finite = state['finite']
    if observed and len(finite) >= 2:
        old = finite[-2]
        _, oldz, oldd = components(old['result'])
        growth = math.log2(n/old['n'])
        if growth > 0:
            if zwidth and oldz:
                zorder = min(12., max(1., (log_fraction(oldz)-log_fraction(zwidth))/growth))
            if remainder and oldd:
                dorder = min(12., max(1., (log_fraction(oldd)-log_fraction(remainder))/growth))
    budget = center*(target-1)/(target+1)
    lz = log_fraction(zwidth)-log_fraction(budget)
    ld = log_fraction(remainder)-log_fraction(budget)
    predicted, steps = 0., 0
    next_n = n
    while next_n < cap:
        next_n = min(2*next_n, cap)
        growth = math.log2(next_n/n)
        a, b = lz-zorder*growth, ld-dorder*growth
        high = max(a, b)
        logsum = high+math.log2(2**(a-high)+2**(b-high)) if math.isfinite(high) else high
        predicted += elapsed*(next_n/n)
        steps += 1
        if logsum <= 0:
            return dict(seconds=predicted, steps=steps, predicted_intervals=next_n,
                        z_order=zorder, remainder_order=dorder, reason='projected_gap')
    return dict(seconds=None, steps=None, z_order=zorder, remainder_order=dorder,
                reason='projected_beyond_mesh_limit')


def run_request(factory, policy, target=F(51, 50), initial=16, cap=4096,
                seconds=30., max_attempts=20, share_nodes=True):
    if policy not in POLICIES:
        raise ValueError('Unknown degree-selection policy')
    if not isinstance(target, F) or target <= 1:
        raise ValueError('Exact rational target greater than one required')
    if type(initial) is not int or type(cap) is not int or not 1 <= initial <= cap <= 16384:
        raise ValueError('Invalid mesh budget')
    if type(max_attempts) is not int or not 1 <= max_attempts <= 64:
        raise ValueError('Invalid attempt budget')
    if isinstance(seconds, bool) or not math.isfinite(seconds) or seconds <= 0:
        raise ValueError('Invalid soft time budget')
    started, cpu = time.perf_counter(), time.process_time()
    problem = factory()
    states, history, predictions = {}, [], []
    final = final_degree = None
    success, stop, committed = False, 'attempt_limit', None

    def state_for(degree):
        if degree not in states:
            v = ShapedVerifier(problem, 'moment') if degree == 2 else CubicVerifier(problem, 'moment')
            if share_nodes and states:
                # Same oracle instance, precision and rational coordinates.
                # Cell caches must remain separate: their derivative orders differ.
                v.node_cache = next(iter(states.values()))['verifier'].node_cache
            states[degree] = dict(degree=degree, verifier=v, n=0, last=None,
                                  last_seconds=0., finite=[])
        return states[degree]

    def choose():
        nonlocal committed
        if policy.startswith('fixed'):
            d = int(policy[5])
            return d if d not in states or states[d]['n'] < cap else None
        if 2 not in states:
            return 2
        if 3 not in states:
            return 3
        available = [d for d in (2, 3) if states[d]['n'] < cap]
        if not available:
            return None
        if policy in ('probe_cubic', 'microprobe_cubic'):
            return 3 if 3 in available else 2
        if policy == 'forecast_guarded':
            # Development exposed a sharp-peak counterexample: on one coarse
            # mesh interval derivative inflation makes nominal orders misleading.
            # Observe actual contraction on two meshes per surviving branch.
            unprobed = [d for d in available if len(states[d]['finite']) < 2]
            if unprobed:
                return min(unprobed, key=lambda d: (states[d]['n'], d))
        if policy == 'forecast_once' and committed in available:
            return committed
        costs = {d: forecast(states[d], target, cap, observed=policy in ('forecast_race', 'forecast_guarded'))
                 for d in available}
        predictions.append(dict(after_attempt=len(history), estimates={str(d): c for d, c in costs.items()}))
        usable = [d for d in available if costs[d]['seconds'] is not None]
        if usable:
            degree = min(usable, key=lambda d: (costs[d]['seconds'], d))
        else:
            # Both forecasts failed. Continue actual certification, alternating
            # toward the less-developed mesh instead of treating failure as proof.
            degree = min(available, key=lambda d: (states[d]['n'], d))
        if policy == 'forecast_once':
            committed = degree
        return degree

    for attempt in range(1, max_attempts+1):
        if time.perf_counter()-started >= seconds:
            stop = 'time_limit'
            break
        degree = choose()
        if degree is None:
            stop = 'mesh_limit'
            break
        state = state_for(degree)
        first = 1 if ((policy == 'microprobe_cubic' and degree == 2) or policy.endswith('_from1')) else initial
        n = first if state['n'] == 0 else min(2*state['n'], cap)
        tick = time.perf_counter()
        row = dict(attempt=attempt, degree=degree, intervals=n)
        try:
            state['verifier'].uniform(n)
            result = state['verifier'].certify()
            passed = bool(gap_met(result, target))
            # Export each attempted certificate in every policy's timed path.
            row.update(status='finite', enclosure=pack(result), goal_met=passed)
            state['last'] = result
            state['finite'].append(dict(n=n, result=result))
            final, final_degree = result, degree
        except ArithmeticError as error:
            state['last'] = None
            row.update(status='arithmetic_unresolved', reason=str(error), goal_met=False)
        state['n'] = n
        state['last_seconds'] = time.perf_counter()-tick
        row['elapsed_seconds'] = state['last_seconds']
        history.append(row)
        if row['goal_met']:
            success, stop = True, 'goal_met'
            break
    # A finite last result remains evidence even if a later attempt is unresolved.
    work = dict(residual_evaluations=sum(s['verifier'].residual_evaluations for s in states.values()),
                derivative_interval_evaluations=sum(s['verifier'].curvature_evaluations for s in states.values()),
                polynomial_cell_visits=sum(s['verifier'].polynomial_cell_visits for s in states.values()),
                problem_forward_calls=getattr(problem, 'forward_calls', None))
    result = dict(policy=policy, target=str(target), success=success, stop_reason=stop,
                  final=pack(final) if final is not None else None, final_degree=final_degree,
                  attempts=len(history), trace=history, predictions=predictions, work=work,
                  wall_ms=(time.perf_counter()-started)*1000,
                  cpu_ms=(time.process_time()-cpu)*1000, shared_residual_nodes=share_nodes)
    return result
