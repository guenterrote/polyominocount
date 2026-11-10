"""
Count fixed polyominoes that are symmetric under reflection in a vertical
axis, with a transfer matrix on the half polyomino.

Reduction: P is invariant under x -> -x (axis through a column of cells,
"odd" case) or x -> 1-x (axis between two columns, "even" case) iff its
right half H is a connected polyomino whose leftmost column is the column
next to / on the axis.  The size of P is

  odd case : |H| + |H minus axis column|   (axis column cells weigh 1, others 2)
  even case: 2 |H|                          (all cells weigh 2)

H is counted with a cell-by-cell transfer matrix in a strip of width W,
rows from top to bottom.  The left edge of the strip is the axis column;
H must touch it.  The frontier holds one connectivity label per column
(0 = empty).  A polyomino is complete when its last component leaves the
frontier; at that moment it is counted (if it touched column 0).

Prototype: exact integer counts, no pruning beyond the size limit.
"""
import sys
from collections import defaultdict


def normalize(f):
    """relabel components in order of first appearance"""
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


def lower_bound(f, touched):
    """cheap lower bound on the number of further cells needed to finish:
    every merge of two components needs a new cell, and reaching column 0
    needs at least (leftmost occupied column) new cells."""
    k = len(set(v for v in f if v))
    lb = max(k - 1, 0)
    if not touched:
        first = next((i for i, v in enumerate(f) if v), 0)
        lb = max(lb, first)
    return lb


STATS = {}


def count_half(nmax, W, w_axis, w_other):
    """counts[n] = number of H (touching column 0, width <= W) with weighted size n"""
    result = defaultdict(int)
    # state: (frontier, touched) -> {weight: count}
    # The empty, not-yet-started state only lives during row 0.
    empty = (0,) * W
    states = {(empty, False): {0: 1}}
    row = 0
    while states:
        for j in range(W):
            wt = w_axis if j == 0 else w_other
            new = defaultdict(lambda: defaultdict(int))
            for (f, touched), gf in states.items():
                started = any(f) or row > 0  # in row 0 the empty state is "not started"
                up = f[j]
                left = f[j - 1] if j > 0 else 0
                # --- option 1: cell (row, j) empty ---
                g = list(f)
                g[j] = 0
                ok = True
                if up and up not in g:
                    # component 'up' leaves the frontier for good
                    if any(g):
                        ok = False  # it can never reconnect: dead
                    else:
                        # polyomino complete
                        if touched:
                            for wgt, c in gf.items():
                                result[wgt] += c
                        ok = False
                if ok and (any(g) or (row == 0)):
                    key = (normalize(g), touched)
                    tgt = new[key]
                    for wgt, c in gf.items():
                        tgt[wgt] += c
                # --- option 2: cell occupied ---
                if not started and row > 0:
                    continue
                g = list(f)
                labels = [l for l in (up, left) if l]
                if labels:
                    lab = labels[0]
                    for other in labels[1:]:
                        if other != lab:
                            g = [lab if v == other else v for v in g]
                else:
                    lab = max(g) + 1
                g[j] = lab
                key = (normalize(g), touched or j == 0)
                tgt = new[key]
                for wgt, c in gf.items():
                    if wgt + wt <= nmax:
                        tgt[wgt + wt] += c
            states = {}
            for k, v in new.items():
                lb = lower_bound(*k) * min(w_axis, w_other)
                v = {wgt: c for wgt, c in v.items() if wgt + lb <= nmax}
                if v:
                    states[k] = v
            STATS['max_states'] = max(STATS.get('max_states', 0), len(states))
        # end of row: the never-started empty state dies
        states.pop((empty, False), None)
        row += 1
    return result


def mirror_counts(nmax):
    W = (nmax + 1) // 2 + 1
    odd = count_half(nmax, W, 1, 2)
    even = count_half(nmax, W, 2, 2)
    return [odd.get(n, 0) + even.get(n, 0) for n in range(nmax + 1)]


if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    c = mirror_counts(nmax)
    for n in range(1, nmax + 1):
        print(n, c[n])
    print("# max states", STATS['max_states'], file=sys.stderr)
