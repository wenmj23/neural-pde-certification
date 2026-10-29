"""Research mesh policies over unchanged quadratic/cubic interval kernels.

Uses cached remainder indicators to choose uniform versus local refinement.
No problem labels, exact error, frequency or reference data select a mesh.
All successful ranges pass the original exact rational gap test.
"""
from fractions import Fraction as F
import math
import heapq
import time

from arb_certification import finite
from arb_cubic_reconstruction import CubicVerifier, NODES, LAGRANGE
from arb_shaped_remainder import ShapedVerifier
from arb_degree_selection import pack, gap_met, components
from arb_reconstruction import point
from flint import arb

POLICIES = ('cubic_uniform16', 'cubic_uniform1', 'cubic_adaptive16',
            'cubic_adaptive1', 'cubic_bulk16', 'quadratic_adaptive16',
            'auto_mesh16', 'auto_mesh1', 'cubic_heap16', 'cubic_heap1',
            'quadratic_adaptive1', 'auto_heap16', 'auto_heap1')


class RecordedCubic(CubicVerifier):
    def prepare(self, nodes):
        cells = super().prepare(nodes)
        self.prepared_nodes, self.prepared_cells = tuple(nodes), cells
        return cells

    def cell(self, a, b):
        """Same arithmetic as CubicVerifier.prepare, for one child interval.

        This additional cell-building path is checked against the original
        kernel on rational meshes before benchmark results are accepted.
        """
        self.guard()
        if not F(0) <= a < b <= F(1): raise ValueError('Invalid child interval')
        if (a, b) in self.cell_cache: return self.cell_cache[(a, b)]
        aa, bb = point(a), point(b); h = finite(bb-aa, 'heap cubic width')
        if h.lower() <= 0: raise ValueError('Unseparated heap interval')
        source = [self._residual(a+(b-a)*t) for t in NODES]
        power = []
        for j in range(4):
            value = arb(0)
            for i in range(4): value += source[i]*LAGRANGE[i][j].numerator/LAGRANGE[i][j].denominator
            power.append(finite(value, 'heap cubic coefficient'))
        k = finite(self.problem.r_fourth(aa.lower().union(bb.upper())), 'heap residual fourth derivative').abs_upper()
        delta = finite(k*h*h*h*h/1944, 'heap cubic remainder').upper()
        self.curvature_evaluations += 1
        cell = dict(a=aa, b=bb, h=h, p=power, delta=delta)
        self.cell_cache[(a, b)] = cell
        return cell

    def heap_to(self, n):
        if type(n) is not int or n < len(self.nodes)-1: raise ValueError('Invalid heap budget')
        self.guard()
        leaves = {(a, b): self.cell(a, b) for a, b in zip(self.nodes[:-1], self.nodes[1:])}
        heap = []
        for (a, b), cell in leaves.items(): heapq.heappush(heap, (-log_score(a, b, cell), a, b))
        while len(leaves) < n:
            _, a, b = heapq.heappop(heap); del leaves[(a, b)]
            m = (a+b)/2
            for left, right in ((a, m), (m, b)):
                cell = self.cell(left, right); leaves[(left, right)] = cell
                heapq.heappush(heap, (-log_score(left, right, cell), left, right))
        self.nodes = sorted({F(0), *[b for _, b in leaves]})


class RecordedQuadratic(ShapedVerifier):
    def prepare(self, nodes):
        cells = super().prepare(nodes)
        self.prepared_nodes, self.prepared_cells = tuple(nodes), cells
        return cells


def log_score(a, b, cell):
    weight = F(1, 4) if a <= F(1, 2) <= b else max(a*(1-a), b*(1-b))
    delta = cell['delta']
    if delta < 0: raise ArithmeticError('Negative mesh indicator')
    value = float(delta.log())+math.log(float((b-a)*weight)) if delta != 0 else -math.inf
    if math.isnan(value) or value == math.inf: raise ArithmeticError('Invalid mesh indicator')
    return value


def indicators(verifier):
    if tuple(verifier.nodes) != verifier.prepared_nodes:
        raise ValueError('Indicators require the last complete mesh')
    logs = []
    for a, b, cell in zip(verifier.nodes[:-1], verifier.nodes[1:], verifier.prepared_cells):
        # Log scaling avoids enormous exact fractions or floating underflow in
        # Gaussian tails. These values select work only, never certify a bound.
        logs.append(log_score(a, b, cell))
    largest = max(logs)
    scores = [math.exp(x-largest) if x != -math.inf else 0. for x in logs] if largest != -math.inf else [0.]*len(logs)
    total = math.fsum(scores); squares = math.fsum(x*x for x in scores)
    concentration = min(1., max(0., total*total/(len(scores)*squares))) if squares else 1.
    return scores, concentration


def bulk_nodes(nodes, scores, cap, mass=F(1, 2)):
    """Bisect a smallest descending set carrying at least half indicator mass.

    The scores are efficiency heuristics; certification is always recomputed.
    """
    if len(scores) != len(nodes)-1 or not 0 < mass <= 1:
        raise ValueError('Invalid bulk-marking input')
    if any(not math.isfinite(x) or x < 0 for x in scores): raise ValueError('Invalid indicator')
    room = cap-(len(nodes)-1)
    if room <= 0: return list(nodes)
    total = math.fsum(scores)
    if not total: return None
    chosen, carried = [], 0.
    for i in sorted(range(len(scores)), key=lambda i: (scores[i], -i), reverse=True):
        chosen.append(i); carried += scores[i]
        if carried >= mass*total or len(chosen) == room: break
    return sorted(list(nodes)+[(nodes[i]+nodes[i+1])/2 for i in chosen])


def next_uniform(n, cap):
    return min(2*n, cap)


def split_all(nodes, cap):
    room = min(cap-(len(nodes)-1), len(nodes)-1)
    chosen = sorted(range(len(nodes)-1), key=lambda i: (nodes[i+1]-nodes[i], -i), reverse=True)[:room]
    return sorted(list(nodes)+[(nodes[i]+nodes[i+1])/2 for i in chosen])


def request(factory, policy, target=F(51, 50), cap=4096, seconds=30., max_attempts=40):
    if policy not in POLICIES: raise ValueError('Unknown mesh policy')
    if not isinstance(target, F) or target <= 1: raise ValueError('Exact gap target > 1 required')
    if type(cap) is not int or not 16 <= cap <= 16384: raise ValueError('Invalid mesh cap')
    if type(max_attempts) is not int or not 1 <= max_attempts <= 64: raise ValueError('Invalid attempts')
    if isinstance(seconds, bool) or not math.isfinite(seconds) or seconds <= 0: raise ValueError('Invalid time cap')
    start, cpu = time.perf_counter(), time.process_time()
    problem = factory()
    degree = 2 if policy.startswith('quadratic') else 3
    verifier = RecordedQuadratic(problem, 'moment') if degree == 2 else RecordedCubic(problem, 'moment')
    n = 1 if policy.endswith('1') and not policy.endswith('16') else 16
    verifier.uniform(n)
    trace = []; final = None; final_nodes = None; success = False; stop = 'attempt_limit'
    action = 'initial_uniform'; locked_mesh = None
    for attempt in range(1, max_attempts+1):
        if time.perf_counter()-start >= seconds:
            stop = 'time_limit'; break
        tick = time.perf_counter()
        row = dict(attempt=attempt, intervals=len(verifier.nodes)-1, action=action)
        diagnostics = None
        try:
            result = verifier.certify()
            passed = bool(gap_met(result, target))
            row.update(status='finite', enclosure=pack(result), goal_met=passed)
            final = result; final_nodes = list(verifier.nodes)
            # Inspect cached cell data only after a failed certificate. Reading
            # diagnostics must not introduce another residual/derivative call.
            if not passed and ('bulk' in policy or (policy.startswith('auto') and locked_mesh != 'uniform')):
                _, zwidth, du = components(result)
                if locked_mesh == 'heap':
                    diagnostics = (None, True, du >= zwidth)
                    row['diagnostics'] = dict(remainder_dominates=bool(du >= zwidth))
                else:
                    scores, concentration = indicators(verifier)
                    localized = concentration < 1/3 and du >= zwidth
                    diagnostics = (scores, localized, du >= zwidth)
                    row['diagnostics'] = dict(effective_cell_fraction=concentration,
                                              remainder_dominates=bool(du >= zwidth), localized=bool(localized))
        except ArithmeticError as error:
            row.update(status='arithmetic_unresolved', reason=str(error), goal_met=False)
        row['elapsed_seconds'] = time.perf_counter()-tick
        trace.append(row)
        if row['goal_met']:
            success, stop = True, 'goal_met'; break
        n = len(verifier.nodes)-1
        if n >= cap:
            stop = 'mesh_limit'; break
        if time.perf_counter()-start >= seconds:
            stop = 'time_limit'; break
        if attempt == max_attempts:
            stop = 'attempt_limit'; break
        # All mesh-construction and marking work remains in the request timer.
        try:
            if 'adaptive' in policy:
                verifier.refine_to(next_uniform(n, cap)); action = 'existing_adaptive'
            elif policy.startswith('cubic_heap'):
                if row['status'] == 'finite' and components(result)[2] < components(result)[1]:
                    verifier.nodes = split_all(verifier.nodes, cap); action = 'heap_global_norm_fallback'
                else:
                    verifier.heap_to(next_uniform(n, cap)); action = 'heap_local'
            elif policy == 'cubic_bulk16':
                next_nodes = bulk_nodes(verifier.nodes, diagnostics[0], cap) if diagnostics and diagnostics[2] else None
                if next_nodes:
                    verifier.nodes = next_nodes; action = 'bulk_half_mass'
                else:
                    verifier.nodes = split_all(verifier.nodes, cap); action = 'global_bisect_norm_or_invalid_fallback'
            elif policy.startswith('auto'):
                if locked_mesh is None and diagnostics:
                    early = 'heap' in policy and n >= 4
                    if n >= 16 or early:
                        threshold = .6 if early and n < 16 else 1/3
                        local = row['diagnostics']['effective_cell_fraction'] < threshold and diagnostics[2]
                        if local or n >= 16:
                            locked_mesh = ('heap' if 'heap' in policy else 'bulk') if local else 'uniform'
                if locked_mesh == 'heap':
                    if diagnostics and diagnostics[2]:
                        verifier.heap_to(next_uniform(n, cap)); action = 'auto_heap_local'
                    else:
                        verifier.nodes = split_all(verifier.nodes, cap); action = 'auto_heap_global_norm_fallback'
                    continue
                if locked_mesh == 'bulk' and diagnostics and diagnostics[2]:
                    next_nodes = bulk_nodes(verifier.nodes, diagnostics[0], cap)
                else: next_nodes = None
                if next_nodes:
                    verifier.nodes = next_nodes; action = 'auto_bulk_half_mass'
                elif locked_mesh == 'bulk':
                    verifier.nodes = split_all(verifier.nodes, cap); action = 'auto_global_bisect'
                else:
                    # Coarse starts 1 -> 4 -> 16; later use ordinary doubling.
                    size = min(4*n, cap) if n < 16 and policy in ('auto_mesh1', 'auto_heap1') else next_uniform(n, cap)
                    verifier.uniform(size); action = 'auto_uniform'
            else:
                verifier.uniform(next_uniform(n, cap)); action = 'uniform'
        except ArithmeticError as error:
            # Do not silently skip an unresolved adaptive step: record it and
            # safely fall back to a new uniform mesh at the next budget.
            verifier.uniform(next_uniform(n, cap)); action = 'adaptive_arithmetic_fallback'
            row['refinement_error'] = str(error)
    work = dict(residual_evaluations=verifier.residual_evaluations,
                derivative_interval_evaluations=verifier.curvature_evaluations,
                polynomial_cell_visits=verifier.polynomial_cell_visits,
                problem_forward_calls=getattr(problem, 'forward_calls', None))
    return dict(policy=policy, target=str(target), degree=degree, success=success, stop_reason=stop,
                final=pack(final) if final is not None else None, attempts=len(trace), trace=trace,
                final_nodes=[str(x) for x in final_nodes] if final_nodes else None,
                locked_mesh=locked_mesh, work=work,
                wall_ms=(time.perf_counter()-start)*1000, cpu_ms=(time.process_time()-cpu)*1000)
