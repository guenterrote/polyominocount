"""
Count fixed polyominoes that are invariant under 180-degree rotation,
with a transfer matrix on the upper half and a "meet in the middle" glue.

The rotation is rho(x, y) = (c - x, d - y).  Four cases, by the centre:
  c=0, d=0 : centre in a cell centre        (middle row 0 exists)
  c=1, d=0 : centre on a vertical edge      (middle row 0 exists)
  c=0, d=1 : centre on a horizontal edge    (no middle row)
  c=1, d=1 : centre on a vertex             (no middle row)
fix(r180) = sum of the four.  (The two edge cases have equal counts by
a 90-degree rotation, but both are computed here as a check.)

U = cells of P with y > d/2 (rows 1, 2, ... above the middle).  Every
component of U must reach row 1.  U is built by an ordinary top-down,
cell-by-cell transfer matrix (non-crossing connectivity labels); its top
row is the first processed row, and after each completed row the frontier
is glued to its rotated copy (and, if d=0, a symmetric middle row R).
Only this glue step sees the twisted connectivity, checked by union-find.

Size: |P| = 2|U| + |R|.
"""
import sys
from collections import defaultdict


def normalize(f):
    m = {}
    out = []
    for v in f:
        if v == 0:
            out.append(0)
        else:
            if v not in m:
                m[v] = len(m) + 1
            out.append(m[v])
    return tuple(out)


class UF:
    def __init__(self):
        self.p = {}

    def find(self, a):
        self.p.setdefault(a, a)
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def union(self, a, b):
        self.p[self.find(a)] = self.find(b)


GLUE_CACHE = {}


def glue(f, xs, c, d, nmax_r):
    key = (f, c, d, nmax_r)
    if key not in GLUE_CACHE:
        GLUE_CACHE[key] = _glue(f, xs, c, d, nmax_r)
    return GLUE_CACHE[key]


def _glue(f, xs, c, d, nmax_r):
    """For frontier f (labels at row 1, columns xs), return {|R|: number of
    middle rows R making P connected}.  If d == 1 there is no middle row."""
    idx = {x: i for i, x in enumerate(xs)}
    labels = set(v for v in f if v)
    res = defaultdict(int)

    def connected(R):
        uf = UF()
        nodes = [("u", a) for a in labels] + [("v", a) for a in labels]
        nodes += [("r", x) for x in R]
        if not nodes:
            return False
        for x in R:
            uf.find(("r", x))
            if x + 1 in R:
                uf.union(("r", x), ("r", x + 1))
            a = f[idx[x]]                  # cell (x, 1) above
            if a:
                uf.union(("r", x), ("u", a))
            b = f[idx[c - x]]              # cell (x, -1) = rho(c - x, 1)
            if b:
                uf.union(("r", x), ("v", b))
        if d == 1:
            # row 1 cell (x,1) is adjacent to (x,0) = rho(c-x, 1)
            for i, x in enumerate(xs):
                a = f[i]
                b = f[idx[c - x]]
                if a and b:
                    uf.union(("u", a), ("v", b))
        root = uf.find(nodes[0])
        return all(uf.find(v) == root for v in nodes)

    if d == 1:
        if labels and connected(set()):
            res[0] += 1
        return res
    # symmetric middle rows: choose cells for x >= c - x, i.e. 2x >= c
    half = [x for x in xs if 2 * x >= c]
    from itertools import combinations
    for k in range(len(half) + 1):
        if (2 * k - (1 if c == 0 and k and 0 in half else 0)) > nmax_r and k > 0:
            # even the cheapest choice (including x=0) is too big
            if 2 * k - 1 > nmax_r:
                break
        for chosen in combinations(half, k):
            R = set()
            for x in chosen:
                R.add(x)
                R.add(c - x)
            if len(R) > nmax_r:
                continue
            if connected(R):
                res[len(R)] += 1
    return res


def count_case(nmax, c, d):
    W = nmax // 2 + 1
    xs = list(range(-W + c, W + 1))     # closed under x -> c - x
    L = len(xs)
    total = defaultdict(int)
    # U empty: P lies in the middle row only (d == 0)
    if d == 0:
        for r, k in glue((0,) * L, xs, c, d, nmax).items():
            total[r] += k
    empty = (0,) * L
    states = {empty: {0: 1}}
    row = 0
    while states:
        for j in range(L):
            new = defaultdict(lambda: defaultdict(int))
            for f, gf in states.items():
                up = f[j]
                left = f[j - 1] if j > 0 else 0
                # empty cell
                g = list(f)
                g[j] = 0
                if not (up and up not in g):     # a component closed off: dead
                    tgt = new[normalize(g)]
                    for u, k in gf.items():
                        tgt[u] += k
                # occupied cell
                if row > 0 and not any(f):
                    continue
                g = list(f)
                labs = [l for l in (up, left) if l]
                if labs:
                    lab = labs[0]
                    if len(labs) == 2 and labs[1] != lab:
                        g = [lab if v == labs[1] else v for v in g]
                else:
                    lab = max(g) + 1
                g[j] = lab
                tgt = new[normalize(g)]
                for u, k in gf.items():
                    if 2 * (u + 1) <= nmax:
                        tgt[u + 1] += k
            states = {k: dict(v) for k, v in new.items() if v}
        states.pop(empty, None)     # U must start in its first row
        row += 1
        # glue: the row just finished is row 1
        for f, gf in states.items():
            gl = glue(f, xs, c, d, nmax - 2 * min(gf))
            for u, k in gf.items():
                for r, m in gl.items():
                    if 2 * u + r <= nmax:
                        total[2 * u + r] += k * m
    return total


def r180_counts(nmax):
    out = [0] * (nmax + 1)
    per_case = {}
    for c in (0, 1):
        for d in (0, 1):
            t = count_case(nmax, c, d)
            per_case[c, d] = t
            for n, k in t.items():
                if 1 <= n <= nmax:
                    out[n] += k
    return out, per_case


if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    out, per = r180_counts(nmax)
    for n in range(1, nmax + 1):
        print(n, out[n])
    for key, t in per.items():
        print("# case c=%d d=%d:" % key, [t.get(n, 0) for n in range(1, nmax + 1)], file=sys.stderr)


# ---------------------------------------------------------------------------
# Version 2: exact width, touching flag, width <= height, pruning.
#
# For each case (c, d) and each half-width X, the strip is xs = [c-X .. X]
# (width w = 2X - c + 1) and P must touch column X (equivalently c - X).
# U has k rows; P has height h = 2k + 1 - d.  Only P with w <= h are
# counted here; the P with w > h of case (c, d) are, after a 90-degree
# rotation, the P with w < h of case (d, c).  Hence
#     fix(r180) = sum over cases of  2 * #(w < h) + #(w == h).

STATS2 = defaultdict(int)


def glue_touch(f, xs, c, d, nmax_r, X, touched):
    """like glue(), but returns {(|R|, R touches column X): count}"""
    key = ("t", f, c, d, nmax_r, X)
    if key not in GLUE_CACHE:
        res = defaultdict(int)
        idx = {x: i for i, x in enumerate(xs)}
        labels = set(v for v in f if v)
        half = [x for x in xs if 2 * x >= c] if d == 0 else []
        from itertools import combinations
        for kk in range(len(half) + 1):
            if 2 * kk - 1 > nmax_r and kk > 0:
                break
            for chosen in combinations(half, kk):
                R = set()
                for x in chosen:
                    R.add(x)
                    R.add(c - x)
                if len(R) > nmax_r:
                    continue
                uf = UF()
                nodes = [("u", a) for a in labels] + [("v", a) for a in labels] + [("r", x) for x in R]
                if not nodes:
                    continue
                for x in R:
                    uf.find(("r", x))
                    if x + 1 in R:
                        uf.union(("r", x), ("r", x + 1))
                    a = f[idx[x]]
                    if a:
                        uf.union(("r", x), ("u", a))
                    b = f[idx[c - x]]
                    if b:
                        uf.union(("r", x), ("v", b))
                if d == 1:
                    for i, x in enumerate(xs):
                        a, b = f[i], f[idx[c - x]]
                        if a and b:
                            uf.union(("u", a), ("v", b))
                root = uf.find(nodes[0])
                if all(uf.find(v) == root for v in nodes):
                    res[len(R), X in R] += 1
        GLUE_CACHE[key] = res
    return GLUE_CACHE[key]


PRUNE_GAPS = False  # tried: valid, but only ~1% fewer states


def gap_bound(f, xs, c):
    """Lower bound on the weight still needed to connect everything.
    Overlay U's frontier (labels u_a at x) with the rotated copy's frontier
    (labels v_a at c - x).  The final P is connected, so every column
    between the outermost occupied ones needs a cell, unless the gap is
    already bridged by old cells: some node u_a or v_a occurring on both
    sides of it.  Each missing column costs at least weight 1."""
    L = len(xs)
    lo, hi = {}, {}
    occ = [False] * L
    for i, x in enumerate(xs):
        j = xs.index(c - x) if False else (L - 1 - i)   # c - xs[i] == xs[L-1-i]
        for node in (("u", f[i]), ("v", f[j])):
            if node[1]:
                occ[i] = True
                if node not in lo:
                    lo[node] = i
                hi[node] = i
    pos = [i for i in range(L) if occ[i]]
    if not pos:
        return 0
    need = 0
    i = pos[0]
    while i <= pos[-1]:
        if occ[i]:
            i += 1
            continue
        g1 = i
        while not occ[i]:
            i += 1
        g2 = i - 1
        if not any(lo[n] < g1 and hi[n] > g2 for n in lo):
            need += g2 - g1 + 1
    return need


def count_case_width(nmax, c, d, X):
    """returns (lt, eq): counts by size of P in case (c,d), half-width X,
    with w < h and w == h respectively"""
    xs = list(range(c - X, X + 1))
    L = len(xs)
    w = L
    lt, eq = defaultdict(int), defaultdict(int)
    edge = (0, L - 1)

    def height(k):
        return 2 * k + 1 - d

    # U empty: P = middle row only; h = 1, so w <= h only if w == 1
    if d == 0 and w == 1:
        eq[1] += 1
    empty = (0,) * L
    states = {(empty, False): {0: 1}}
    row = 0
    kmin = 0
    while height(kmin) < w:
        kmin += 1
    while states:
        for j in range(L):
            new = defaultdict(lambda: defaultdict(int))
            for (f, t), gf in states.items():
                up = f[j]
                left = f[j - 1] if j > 0 else 0
                g = list(f)
                g[j] = 0
                if not (up and up not in g):
                    tgt = new[(normalize(g), t)]
                    for u, k in gf.items():
                        tgt[u] += k
                if row > 0 and not any(f):
                    continue
                g = list(f)
                labs = [l for l in (up, left) if l]
                if labs:
                    lab = labs[0]
                    if len(labs) == 2 and labs[1] != lab:
                        g = [lab if v == labs[1] else v for v in g]
                else:
                    lab = max(g) + 1
                g[j] = lab
                tgt = new[(normalize(g), t or j in edge)]
                for u, k in gf.items():
                    if 2 * (u + 1) <= nmax:
                        tgt[u + 1] += k
            # pruning: rows still needed and distance to the edge columns
            states = {}
            rows_done = row + 1  # the current row may already hold its cell
            for (f, t), gf in new.items():
                need_rows = max(0, kmin - rows_done)
                lb = 2 * need_rows
                if PRUNE_GAPS:
                    lb = max(lb, gap_bound(f, xs, c))
                if not t:
                    occ = [i for i, v in enumerate(f) if v]
                    if occ:
                        dist = min(min(i, L - 1 - i) for i in occ)
                        lb = max(lb, 2 * dist)
                gf2 = {u: k for u, k in gf.items() if 2 * u + lb <= nmax}
                if gf2:
                    states[(f, t)] = gf2
            STATS2['states'] = max(STATS2['states'], len(states))
        states.pop((empty, False), None)
        row += 1
        h = height(row)
        if h < w:
            continue
        target = eq if h == w else lt
        for (f, t), gf in states.items():
            umin = min(gf)
            gl = glue_touch(f, xs, c, d, nmax - 2 * umin, X, t)
            for u, k in gf.items():
                for (r, rt), m in gl.items():
                    if (t or rt) and 2 * u + r <= nmax:
                        target[2 * u + r] += k * m
    return lt, eq


def r180_counts_v2(nmax):
    out = [0] * (nmax + 1)
    for c in (0, 1):
        for d in (0, 1):
            for X in range(c, nmax // 2 + 2):
                if 2 * X - c + 1 > nmax:
                    break
                lt, eq = count_case_width(nmax, c, d, X)
                for n, k in lt.items():
                    if n <= nmax:
                        out[n] += 2 * k
                for n, k in eq.items():
                    if n <= nmax:
                        out[n] += k
    return out
