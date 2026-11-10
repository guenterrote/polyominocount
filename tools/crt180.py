"""Combine diag180 results modulo several moduli (file results180.py) by CRT.

usage: python3 tools/crt180.py [results180.py] [bits]
Prints n, R180C, R180M, R180V, total.  A value is accepted only if its CRT
reconstruction is below 2**bits (default 72); a wrong run would give a
value of the size of the product of the moduli.
"""
import sys
from math import prod

fn = sys.argv[1] if len(sys.argv) > 1 else "results180.py"
bits = int(sys.argv[2]) if len(sys.argv) > 2 else 72
env = {}
exec(open(fn).read(), env)
R = env["Results"]
ns = sorted({n for n, m, t in R})
mods = sorted({m for n, m, t in R})
M = prod(mods)
for n in ns:
    vals = []
    for t in ("R180C", "R180M", "R180V"):
        x = 0
        for m in mods:
            r = R.get((n, m, t))
            if r is None:
                sys.exit(f"missing Results[{n}, {m}, {t!r}]")
            Mi = M // m
            x += r * Mi * pow(Mi, -1, m)
        x %= M
        if x >= 2 ** bits:
            sys.exit(f"n={n} {t}: CRT value has {x.bit_length()} bits > {bits}: inconsistent runs?")
        vals.append(x)
    print(n, *vals, sum(vals))
