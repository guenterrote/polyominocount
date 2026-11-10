"""
FIXED polyominoes with the same diagonal-row transfer matrix and pruning as
rot180_diag.py, only to compare state counts with the 180-degree version.

Rows u = x + y = 0, 1, 2, ... (P starts in row 0), frontier over v in
[0, L-1], P touches v = 0 and v = L-1.  Cells have v = u + a (mod 2), a in
{0, 1} (two translation classes of the normalised position).
Counting only vspan <= uspan: fixed = sum 2 #(vspan < uspan) + #(=).
"""
import sys
from collections import defaultdict

STATS = defaultdict(int)


def canon(f, parent):
    def find(a):
        while parent.get(a, a) != a:
            a = parent[a]
        return a
    m, out = {}, []
    for v in f:
        if not v:
            out.append(0)
            continue
        r = find(v)
        if r not in m:
            m[r] = len(m) + 1
        out.append(m[r])
    return tuple(out)


def bound(f, tl, tr):
    occ = [i for i, v in enumerate(f) if v]
    if not occ:
        return 0
    parent = {}

    def find(a):
        while parent.get(a, a) != a:
            a = parent[a]
        return a
    ncls = len(set(f[i] for i in occ))
    cost = 0
    for wgt, x, y in sorted((b - a - 1, f[a], f[b]) for a, b in zip(occ, occ[1:])):
        if ncls == 1:
            break
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry
            cost += wgt
            ncls -= 1
    if not tl:
        cost += occ[0]
    if not tr:
        cost += len(f) - 1 - occ[-1]
    return cost


def count_width(nmax, L, a):
    lt, eq = defaultdict(int), defaultdict(int)
    if L == 1:
        if a == 0:
            eq[1] += 1
        return lt, eq
    kmin = L            # rows needed: uspan >= vspan
    states = {((0,) * L, False, False): {0: 1}}
    u = 0
    while states:
        cells = [i for i in range(L) if (i - u - a) % 2 == 0]
        cur = states
        for ci, i in enumerate(cells):
            last = ci == len(cells) - 1
            nxt = defaultdict(lambda: defaultdict(int))
            for (f, tl, tr), gf in cur.items():
                empty = not any(f)
                for occ in (False, True):
                    if occ and u > 0 and empty:
                        continue
                    f2 = list(f)
                    parent = {}

                    def find(x):
                        while parent.get(x, x) != x:
                            x = parent[x]
                        return x
                    if occ:
                        new = max(f2) + 1
                        f2[i] = new
                        for nb in (i - 1, i + 1):
                            if 0 <= nb < L and f2[nb]:
                                rx, ry = find(f2[nb]), find(new)
                                if rx != ry:
                                    parent[rx] = ry
                    status = 'ok'
                    for j in ([i - 1] if i >= 1 else []) + ([i + 1] if last and i + 1 < L else []):
                        if status != 'ok' or not f2[j]:
                            continue
                        r = find(f2[j])
                        f2[j] = 0
                        if not any(v and find(v) == r for v in f2):
                            status = 'complete' if not any(f2) else 'dead'
                    if status == 'dead':
                        continue
                    ntl = tl or (occ and i == 0)
                    ntr = tr or (occ and i == L - 1)
                    if status == 'complete':
                        if tl and tr:
                            h = u            # rows 0 .. u-1
                            if h >= L:
                                tgt = eq if h == L else lt
                                for wt, cnt in gf.items():
                                    tgt[wt] += cnt
                        continue
                    nf = canon(f2, parent)
                    add = 1 if occ else 0
                    tg = nxt[(nf, ntl, ntr)]
                    for wt, cnt in gf.items():
                        if wt + add <= nmax:
                            tg[wt + add] += cnt
            cur = {}
            for (f, tl, tr), gf in nxt.items():
                if not any(f):
                    if u == 0 and not last:
                        cur[(f, tl, tr)] = dict(gf)
                    continue
                lb = max(0, kmin - u - 1)
                lb = max(lb, bound(f, tl, tr))
                g2 = {w: c for w, c in gf.items() if w + lb <= nmax}
                if g2:
                    cur[(f, tl, tr)] = g2
            STATS['states'] = max(STATS['states'], len(cur))
        states = cur
        u += 1
    return lt, eq


def fixed_counts(nmax):
    out = [0] * (nmax + 1)
    for L in range(1, nmax + 1):
        for a in (0, 1):
            lt, eq = count_width(nmax, L, a)
            for n, k in lt.items():
                out[n] += 2 * k
            for n, k in eq.items():
                out[n] += k
    return out


if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    out = fixed_counts(nmax)
    for n in range(1, nmax + 1):
        print(n, out[n])
    print("# max states", STATS['states'], file=sys.stderr)
