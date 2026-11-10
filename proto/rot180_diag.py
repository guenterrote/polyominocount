"""
180-degree symmetric polyominoes: outward transfer matrix along DIAGONAL
rows (in the spirit of Barequet and Ben-Shachar).

Coordinates u = x + y, v = x - y (cells have u = v mod 2).  Neighbours of
(u, v) are (u +- 1, v +- 1): a diagonal row u = const has no internal
adjacency, and each cell touches at most two cells of the previous row.

rho(x, y) = (c - x, d - y) becomes (u, v) -> (s - u, t - v) with s = c + d,
t = c - d.  If s is even, the diagonal row u = s/2 is mapped to itself
(middle row R, v -> t - v).  If s is odd, the first row u1 = (s+1)/2 of the
half H is adjacent to its own image (row u1 - 1): twin links inside row u1.

The frontier is an array over v in [vlo, vhi]; positions of parity u hold
row u (once processed), the others hold row u - 1 until they are retired.
State bookkeeping (classes, partners, ghosts) as in rot180_out.py.

Width: P spans v in [t - vhi, vhi] and touches vhi; vspan = 2 vhi - t + 1.
Height: uspan = 2k + 1 - (s mod 2) with k rows of H.  Counting only
vspan <= uspan: fix(r180) = sum over cases of 2 #(vspan < uspan) + #(=)
(the reflection (x, y) -> (x, -y) swaps u and v and preserves each case).
"""
import sys
from collections import defaultdict
from itertools import combinations
from rot180_out import State, completion_possible

STATS = defaultdict(int)


def retire(st, i):
    """remove the frontier entry at index i; returns 'ok'/'dead'/'complete'"""
    lab = st.f[i]
    if not lab:
        return 'ok'
    st.f[i] = 0
    r = st.find(lab)
    if any(v and st.find(v) == r for v in st.f):
        return 'ok'
    p = st.rho(r)
    if p == r:
        return 'complete' if not any(st.f) else 'dead'
    if isinstance(p, tuple):
        return 'dead'
    st.partner[st.find(p)] = 'g'
    return 'ok'


def gap_bound(f, p):
    """new cells needed to join different super-classes: a path between
    frontier positions v1 < v2 changes v by +-1 per cell"""
    part = dict(p)
    sup = {lab: (min(lab, q) if isinstance(q, int) else lab) for lab, q in part.items()}
    occ = [i for i, v in enumerate(f) if v]
    if not occ:
        return 0
    lo, hi = {}, {}
    for i in occ:
        s_ = sup[f[i]]
        lo.setdefault(s_, i)
        hi[s_] = i
    if len(lo) == 1:
        return 0
    need = 0
    for a, b in zip(occ, occ[1:]):
        if b - a > 1 and not any(lo[s_] <= a and hi[s_] >= b for s_ in lo):
            need += b - a - 1
    return need


PRUNE_EXACT = True
REFLECT = True


def exact_bound(f, p, touched, kink=None):
    """Cells still needed (Rote's observation: on the diagonal lattice the
    connections do not overlap).  A path between frontier positions i < j
    needs one cell for every v strictly between them, and paths for
    different gaps use disjoint v ranges, so the cheapest way to join the
    super-classes is a minimum spanning tree over the gaps between
    consecutive occupied positions, with weight (j - i - 1).  Reaching an
    edge of the strip (if not yet touched) costs the distance from the
    outermost occupied position, in v columns no gap uses."""
    part = dict(p)
    sup = {lab: (min(lab, q) if isinstance(q, int) else lab) for lab, q in part.items()}
    occ = [i for i, v in enumerate(f) if v]
    if not occ:
        return 0
    parent = {}

    def find(a):
        while parent.get(a, a) != a:
            a = parent[a]
        return a
    edges = sorted((b - a - 1 + (b == kink), sup[f[a]], sup[f[b]]) for a, b in zip(occ, occ[1:]))
    ncls = len(set(sup[f[i]] for i in occ))
    cost = 0
    for wgt, x, y in edges:
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry
            cost += wgt
            ncls -= 1
            if ncls == 1:
                break
    if not touched:
        cost += min(occ[0], len(f) - 1 - occ[-1])
    return cost


SELF_BOUND = True
KINK = True


def self_extra(f, p, kink=None):
    """Extra cells needed to create the self class when there is none.
    Present-partner classes nest as A1 < ... < Ak < rho Ak < ... < rho A1.
    With A_i ~ rho A_i, the top-level gaps form cycles: ring i (left path
    A_i..A_{i+1} plus right path rho A_{i+1}..rho A_i) and the middle
    path A_k..rho A_k.  The MST omits the largest gap of every cycle; but
    some X must get joined to rho X, which needs one of these cycles to be
    closed completely, i.e. its largest gap paid as well."""
    part = dict(p)
    if any(q == 's' for q in part.values()):
        return 0
    runs = []
    for i, v in enumerate(f):
        if v and (i == 0 or not f[i - 1]):
            runs.append([i, i, v])
        elif v:
            runs[-1][1] = i
    lastrun = {}
    for j, r in enumerate(runs):
        lastrun[r[2]] = j
    top, gaps = [], []
    j = 0
    while j < len(runs):
        lab = runs[j][2]
        top.append(lab)
        jl = lastrun[lab]
        if jl + 1 < len(runs):
            gaps.append(runs[jl + 1][0] - runs[jl][1] - 1 + (runs[jl + 1][0] == kink))
        j = jl + 1
    cross = [t for t, lab in enumerate(top) if part[lab] != 'g']
    m = len(cross)
    if m < 2:
        return 0
    k = m // 2
    best = None
    for i in range(k):
        if i < k - 1:
            seg = gaps[cross[i]:cross[i + 1]] + gaps[cross[m - 2 - i]:cross[m - 1 - i]]
        else:
            seg = gaps[cross[k - 1]:cross[k]]
        mx = max(seg)
        best = mx if best is None else min(best, mx)
    return best


def count_case_width(nmax, c, d, vhi):
    s, t = c + d, c - d
    vlo = t - vhi
    L = vhi - vlo + 1
    vspan = L
    lt, eq = defaultdict(int), defaultdict(int)
    e = s % 2
    if L == 1:          # a single diagonal: only the monomino (centre in a cell)
        if (c, d) == (0, 0):
            eq[1] += 1
        return lt, eq

    def uspan(k):
        return 2 * k + 1 - e

    kmin = 0
    while uspan(kmin) < vspan:
        kmin += 1

    states = defaultdict(lambda: defaultdict(int))
    if e == 0:
        um = s // 2
        pos = [i for i in range(L) if (vlo + i - um) % 2 == 0]
        half = [i for i in pos if 2 * (vlo + i) >= t]
        for kk in range(len(half) + 1):
            if 2 * kk - 1 > nmax:
                break
            for chosen in combinations(half, kk):
                cells = set()
                for i in chosen:
                    cells.add(i)
                    cells.add(L - 1 - i)         # v -> t - v
                if not cells or len(cells) > nmax:
                    continue
                f = [0] * L
                partner = {}
                lab = {}
                for i in sorted(cells):
                    lab[i] = len(lab) + 1
                    f[i] = lab[i]
                for i in cells:
                    j = L - 1 - i
                    partner[lab[i]] = 's' if j == i else lab[j]
                st = State(f, partner)
                key = st.canon() + ((0 in cells or (L - 1) in cells),)
                states[key][len(cells)] += 1
        u0 = um + 1
    else:
        states[((0,) * L, (), False)][0] = 1
        u0 = (s + 1) // 2

    def record(k_top, gf, touched):
        h = uspan(k_top)
        if not touched or h < vspan:
            return
        tgt = eq if h == vspan else lt
        for wt, cnt in gf.items():
            tgt[wt] += cnt

    k = 1
    while states:
        u = u0 + k - 1
        cells_u = [i for i in range(L) if (vlo + i - u) % 2 == 0]
        cur = dict(states)
        first_twin = (e == 1 and k == 1)
        for ci, i in enumerate(cells_u):
            nxt = defaultdict(lambda: defaultdict(int))
            for (f, p, tch), gf in cur.items():
                for occ in (False, True):
                    if occ and k > 1 and not any(f):
                        continue
                    st = State(f, dict(p))
                    status = 'ok'
                    if occ:
                        new = max([v for v in st.f if v] + [0] + list(st.partner)) + 1
                        st.partner[new] = 'g'
                        old = st.f[i]
                        st.f[i] = new
                        for nb in (i - 1, i + 1):
                            if 0 <= nb < L and st.f[nb] and nb != i:
                                st.union(new, st.f[nb])
                        if first_twin:
                            # (u1, v) ~ (u1 - 1, v +- 1) = rho(u1, t - v -+ 1)
                            for jj in (L - 1 - i - 1, L - 1 - i + 1):
                                if jj == i:
                                    st.union(new, st.rho(new))
                                elif 0 <= jj < L and jj < i and st.f[jj]:
                                    st.union(new, st.rho(st.f[jj]))
                    # retire the row-(u-1) entry at i - 1
                    if i - 1 >= 0:
                        status = retire(st, i - 1)
                    if status == 'ok' and ci == len(cells_u) - 1 and i + 1 < L:
                        status = retire(st, i + 1)
                    if status == 'dead':
                        continue
                    if status == 'complete':
                        record(k - 1, gf, tch)
                        continue
                    nf, np_ = st.canon()
                    nt = tch or (occ and (i == 0 or i == L - 1))
                    add = 2 if occ else 0
                    tgt = nxt[(nf, np_, nt)]
                    for wt, cnt in gf.items():
                        if wt + add <= nmax:
                            tgt[wt + add] += cnt
            cur = {}
            for (f, p, tch), gf in nxt.items():
                if not any(f):
                    if first_twin and ci < len(cells_u) - 1:
                        cur[(f, p, tch)] = dict(gf)   # first row not started yet
                    continue
                lb = 2 * max(0, kmin - k)
                if not tch:
                    occ_i = [j for j, v in enumerate(f) if v]
                    dist = min(min(j, L - 1 - j) for j in occ_i)
                    lb = max(lb, 2 * dist)
                if not first_twin or ci == len(cells_u) - 1:
                    if not completion_possible(p):
                        continue
                    if PRUNE_EXACT:
                        # kink: position i decided empty in row u, i+1 holds a
                        # row u-1 cell b; a path from the left reaches b only
                        # via (u, i+2), crossing column i+1 above b: +1 cell
                        kk = (i + 1 if KINK and ci < len(cells_u) - 1 and not f[i]
                              and i + 1 < L and f[i + 1] else None)
                        eb = exact_bound(f, p, tch, kk)
                        if SELF_BOUND:
                            eb += self_extra(f, p, kk)
                        lb = max(lb, 2 * eb)
                    else:
                        lb = max(lb, 2 * gap_bound(f, p))
                gf2 = {wt: cnt for wt, cnt in gf.items() if wt + lb <= nmax}
                if gf2:
                    cur[(f, p, tch)] = gf2
            STATS['states'] = max(STATS['states'], len(cur))
        if REFLECT and t == 0:
            # Between rows the frontier holds one full diagonal; for t = 0 the
            # reflection v -> -v (x <-> y) commutes with rho and maps the strip
            # to itself, so a state and its mirror image have the same future.
            merged = defaultdict(lambda: defaultdict(int))
            for (f, p, tch), gf in cur.items():
                key = (f, p, tch)
                st = State(tuple(reversed(f)), dict(p))
                rf, rp = st.canon()
                rkey = (rf, rp, tch)
                key = min(key, rkey, key=repr)
                for wt, cnt in gf.items():
                    merged[key][wt] += cnt
            cur = {k: dict(v) for k, v in merged.items()}
        states = cur
        k += 1
    return lt, eq


def r180_counts(nmax):
    out = [0] * (nmax + 1)
    for c in (0, 1):
        for d in (0, 1):
            t = c - d
            vhi = -(-t // 2) if t >= 0 else (t + 1) // 2
            while 2 * vhi - t + 1 <= 2 * nmax:
                if 2 * vhi - t + 1 >= 1:
                    lt, eq = count_case_width(nmax, c, d, vhi)
                    for n, k in lt.items():
                        if 1 <= n <= nmax:
                            out[n] += 2 * k
                    for n, k in eq.items():
                        if 1 <= n <= nmax:
                            out[n] += k
                vhi += 1
    return out


if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    out = r180_counts(nmax)
    for n in range(1, nmax + 1):
        print(n, out[n])
    print("# max states", STATS['states'], file=sys.stderr)
