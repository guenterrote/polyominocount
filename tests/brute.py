"""Independent brute-force reference: grow free polyominoes size by size,
canonicalise under translation (fixed) and under D4 (free)."""
import sys

def norm(cells):
    mx = min(x for x, y in cells); my = min(y for x, y in cells)
    return frozenset((x - mx, y - my) for x, y in cells)

def d4(cells):
    ts = [lambda x, y: (x, y), lambda x, y: (-y, x), lambda x, y: (-x, -y), lambda x, y: (y, -x),
          lambda x, y: (-x, y), lambda x, y: (x, -y), lambda x, y: (y, x), lambda x, y: (-y, -x)]
    return [norm([t(x, y) for x, y in cells]) for t in ts]

def main(N):
    fixed = {norm([(0, 0)])}
    for n in range(1, N + 1):
        free = {min(d4(p), key=lambda s: sorted(s)) for p in fixed}
        print(n, len(fixed), len(free))
        nxt = set()
        for p in fixed:
            for x, y in p:
                for c in ((x+1, y), (x-1, y), (x, y+1), (x, y-1)):
                    if c not in p:
                        nxt.add(norm(p | {c}))
        fixed = nxt

main(int(sys.argv[1]))
