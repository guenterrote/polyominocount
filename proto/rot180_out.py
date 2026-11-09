"""
180-degree symmetric polyominoes: transfer matrix running OUTWARD from the
middle, on the quotient P / rho.

rho(x, y) = (c - x, d - y), cases (c, d) in {0,1}^2 as in rot180_tm.py.
We build P_k = R  u  U_k  u  rho(U_k), where R is the (symmetric) middle row
(only if d = 0) and U_k = cells of U in rows 1..k.  The frontier is row k of
U.  Components of P_k are permuted by rho.  The state records, for the
components that meet the frontier ("present" classes):
  * the frontier labels,
  * for each present class its partner rho(class): itself ("s"), another
    present class, or an absent class ("g", a ghost that has no frontier
    cell; it can still merge, as the mirror image of merges of its partner).
Every union(a, b) is accompanied by union(rho a, rho b).

P is complete when the frontier empties; then it must consist of one
self-partnered class.  Size: |P| = |R| + 2 |U|.

Exact width (strip xs = [c-X .. X], P touches column X) and width <= height
as in rot180_tm.py: fix(r180) = sum over cases of 2 #(w < h) + #(w == h).
"""
import sys
from collections import defaultdict
from itertools import combinations


class State:
    """mutable working copy: frontier labels (0 = empty) and a partner map
    over nodes; nodes are positive ints (classes), ghosts are ('g', i)."""

    def __init__(self, f, partner):
        self.f = list(f)
        self.parent = {}
        self.partner = dict(partner)      # class -> class / 's' / 'g'

    def find(self, a):
        while self.parent.get(a, a) != a:
            a = self.parent[a]
        return a

    def rho(self, a):
        a = self.find(a)
        if isinstance(a, tuple):              # ghost of i
            return self.find(a[1])
        p = self.partner[a]
        if p == 's':
            return a
        if p == 'g':
            return self.find(('g', a))
        return self.find(p)

    def _link(self, a, b):
        a, b = self.find(a), self.find(b)
        if a == b:
            return
        # prefer a real class as representative
        if isinstance(a, tuple):
            a, b = b, a
        self.parent[b] = a

    def union(self, a, b):
        ra, rb = self.rho(a), self.rho(b)
        self._link(a, b)
        self._link(ra, rb)

    def canon(self):
        """normalize labels by first appearance; recompute partners"""
        m = {}
        out = []
        for v in self.f:
            if v == 0:
                out.append(0)
                continue
            r = self.find(v)
            if r not in m:
                m[r] = len(m) + 1
            out.append(m[r])
        partner = {}
        for r, i in m.items():
            q = self.rho(r)
            if q == r:
                partner[i] = 's'
            elif q in m:
                partner[i] = m[q]
            else:
                partner[i] = 'g'
        return tuple(out), tuple(sorted(partner.items()))


def step_cell(f, partner, j, occupied, L, twin_rule, row_cells):
    """process cell j of the current row.  Returns (status, newf, newpartner)
    status: 'ok', 'dead', 'complete'.  row_cells: dict index->label of cells
    of the current row already decided (used for the d=1 twin rule)."""
    st = State(f, dict(partner))
    up = f[j]
    left = f[j - 1] if j > 0 else 0
    if not occupied:
        st.f[j] = 0
        if up:
            r = st.find(up)
            still = any(v and st.find(v) == r for v in st.f)
            if not still:
                p = st.rho(r)
                if p == r:
                    # a self-partnered class is closed for good
                    if any(st.f):
                        return 'dead', None, None
                    return 'complete', None, None
                if isinstance(p, tuple):
                    return 'dead', None, None      # class and its image both closed
                # partner p is present: r becomes the ghost of p
                st.partner[st.find(p)] = 'g'
        nf, np_ = st.canon()
        return 'ok', nf, np_
    # occupied
    new = max([v for v in st.f if v] + [0] + [k for k in st.partner]) + 1
    st.partner[new] = 'g'
    st.f[j] = new
    if up:
        st.union(new, up)
    if left:
        st.union(new, left)
    if twin_rule is not None:
        jj = twin_rule(j)
        if jj == j:
            st.union(new, st.rho(new))
        elif jj in row_cells:
            st.union(new, st.rho(st.f[jj]))
    nf, np_ = st.canon()
    return 'ok', nf, np_


def count_case_width(nmax, c, d, X):
    xs = list(range(c - X, X + 1))
    L = len(xs)
    w = L
    edge_cols = (0, L - 1)
    lt, eq = defaultdict(int), defaultdict(int)

    def height(k):
        return 2 * k + 1 - d

    # ---- initial states: the middle row R (d = 0) or nothing (d = 1)
    states = defaultdict(lambda: defaultdict(int))
    if d == 0:
        half = [i for i in range(L) if 2 * xs[i] >= c]
        for k in range(len(half) + 1):
            if 2 * k - 1 > nmax:
                break
            for chosen in combinations(half, k):
                cells = set()
                for i in chosen:
                    cells.add(i)
                    cells.add(L - 1 - i)
                if not cells or len(cells) > nmax:
                    continue
                # runs of R; rho maps run [a,b] to [L-1-b, L-1-a]
                f = [0] * L
                lab = 0
                runs = {}
                for i in range(L):
                    if i in cells:
                        if i == 0 or (i - 1) not in cells:
                            lab += 1
                            start = i
                        f[i] = lab
                        runs[lab] = (start, i)
                start_of = {v: k_ for k_, v in runs.items()}
                partner = {}
                for lb, (a, b) in runs.items():
                    img = (L - 1 - b, L - 1 - a)
                    q = start_of[img]
                    partner[lb] = 's' if q == lb else q
                st = State(f, partner)
                key = st.canon() + ((0 in cells or (L - 1) in cells),)
                states[key][len(cells)] += 1
    else:
        states[((0,) * L, (), False)][0] = 1

    kmin = 0
    while height(kmin) < w:
        kmin += 1
    twin = (lambda j: L - 1 - j) if d == 1 else None

    def record(k_top, gf, touched):
        h = height(k_top)
        if not touched or h < w:
            return
        tgt = eq if h == w else lt
        for wt, cnt in gf.items():
            tgt[wt] += cnt

    row = 1
    while states:
        # process row `row` cell by cell; row_cells tracks decided cells
        cur = {(f, p, t, ()): gf for (f, p, t), gf in states.items()}
        for j in range(L):
            nxt = defaultdict(lambda: defaultdict(int))
            for (f, p, t, rc), gf in cur.items():
                rcd = dict.fromkeys(rc)
                for occ in (False, True):
                    if occ and row > 1 and not any(f):
                        continue
                    if occ and d == 1 and row == 1 and False:
                        pass
                    rc_now = dict(rcd)
                    status, nf, np_ = step_cell(f, dict(p), j, occ, L,
                                                twin if row == 1 else None, rc_now)
                    if status == 'dead':
                        continue
                    if status == 'complete':
                        record(row - 1, gf, t)
                        continue
                    nt = t or (occ and j in edge_cols)
                    nrc = rc + (j,) if occ else rc
                    tgt = nxt[(nf, np_, nt, nrc if (d == 1 and row == 1) else ())]
                    add = 2 if occ else 0
                    for wt, cnt in gf.items():
                        if wt + add <= nmax:
                            tgt[wt + add] += cnt
            # pruning: rows still needed for h >= w, distance to the edge
            cur = {}
            for key, gf in nxt.items():
                f, p, t, rc = key
                lb = 2 * max(0, kmin - row)
                if PRUNE_MERGE:
                    lb = max(lb, 2 * merge_bound(f, p))
                if not t:
                    occ_i = [i for i, v in enumerate(f) if v]
                    if occ_i:
                        dist = min(min(i, L - 1 - i) for i in occ_i)
                        lb = max(lb, 2 * dist)
                gf2 = {wt: cnt for wt, cnt in gf.items() if wt + lb <= nmax}
                if gf2:
                    cur[key] = gf2
            STATS['states'] = max(STATS['states'], len(cur))
        # end of row: an all-empty frontier that never started (d=1, row 1)
        states = {}
        for (f, p, t, rc), gf in cur.items():
            if not any(f):
                continue          # (complete cases were recorded already)
            states[(f, p, t)] = gf
        row += 1
    return lt, eq


STATS = defaultdict(int)
PRUNE_MERGE = True


def merge_bound(f, p):
    """Number of new cells needed so that all frontier classes can end in
    one component.  A class and its present partner form a super-class.
    Different super-classes can only be joined by paths above the
    frontier, so each empty frontier column in a gap that separates
    super-classes (no super-class on both sides) needs a new cell."""
    part = dict(p)
    sup = {}
    for lab, q in part.items():
        sup[lab] = min(lab, q) if isinstance(q, int) else lab
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


def r180_counts(nmax):
    out = [0] * (nmax + 1)
    for c in (0, 1):
        for d in (0, 1):
            for X in range(c, nmax // 2 + 2):
                if 2 * X - c + 1 > nmax:
                    break
                lt, eq = count_case_width(nmax, c, d, X)
                for n, k in lt.items():
                    if 1 <= n <= nmax:
                        out[n] += 2 * k
                for n, k in eq.items():
                    if 1 <= n <= nmax:
                        out[n] += k
    return out


if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    out = r180_counts(nmax)
    for n in range(1, nmax + 1):
        print(n, out[n])
    print("# max states", STATS['states'], file=sys.stderr)
