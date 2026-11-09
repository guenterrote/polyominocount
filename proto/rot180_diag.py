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
                    lb = max(lb, 2 * gap_bound(f, p))
                gf2 = {wt: cnt for wt, cnt in gf.items() if wt + lb <= nmax}
                if gf2:
                    cur[(f, p, tch)] = gf2
            STATS['states'] = max(STATS['states'], len(cur))
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
