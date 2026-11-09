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


def glue(f, xs, c, d, nmax_r):
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
    for mask in range(1 << len(half)):
        R = set()
        for k, x in enumerate(half):
            if mask >> k & 1:
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
            gl = glue(f, xs, c, d, nmax)
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
