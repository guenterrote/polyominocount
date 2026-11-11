"""Combine diag180 results modulo several moduli (file results180.py) by CRT.

usage: python3 tools/crt180.py [results180.py] [bits] [-b PREFIX]

Prints
  1. a table n, R180C, R180M, R180V, total;
  2. for each column the bare comma-separated list of numbers, as wanted for
     the OEIS data field (with the OEIS number where the column is an OEIS
     sequence), and the same list restricted to the n for which the column
     can be nonzero (R180C: odd n; R180M, R180V: even n);
  3. with -b PREFIX, b-files PREFIX-R180C.txt, ... in OEIS format
     (lines "n a(n)").

A value is accepted only if its CRT reconstruction is below 2**bits
(default 72); a wrong run would give a value of the size of the product of
the moduli.  A result file with modulus 2**128 (128-bit build) needs no
second modulus.

The columns count FIXED polyominoes that are invariant under the rotation by
180 degrees (whatever other symmetries they have), by the type of centre:
    fix(r180) = 4*A006747 + 2*A056877 + 2*A056878 + 2*A144553 + A142886.
The OEIS sequences A351615 / A234008 / A351616 (centre in a cell / on an edge
/ at a vertex) count FREE polyominoes with exactly this symmetry, so they
are not our columns; they follow once the other symmetry classes are known.
"""
import sys
from math import prod

# column -> (description, OEIS number or None, parity of n with nonzero terms)
COLUMNS = {
    "R180C": ("fixed polyominoes invariant under r180, centre in a cell", None, 1),
    "R180M": ("fixed polyominoes invariant under r180, centre on an edge midpoint", None, 0),
    "R180V": ("fixed polyominoes invariant under r180, centre at a vertex", None, 0),
    "total": ("fixed polyominoes invariant under r180", None, None),
}


def combine(fn, bits):
    env = {}
    exec(open(fn).read(), env)
    R = env["Results"]
    ns = sorted({n for n, m, t in R})
    mods = sorted({m for n, m, t in R})
    M = prod(mods)
    table = {}
    for n in ns:
        vals = {}
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
            vals[t] = x
        vals["total"] = sum(vals.values())
        table[n] = vals
    return ns, mods, table


def main():
    args = sys.argv[1:]
    prefix = None
    if "-b" in args:
        i = args.index("-b")
        prefix = args[i + 1]
        del args[i:i + 2]
    fn = args[0] if len(args) > 0 else "results180.py"
    bits = int(args[1]) if len(args) > 1 else 72
    ns, mods, table = combine(fn, bits)

    print("# moduli:", ", ".join(str(m) for m in mods))
    print("# n R180C R180M R180V total")
    for n in ns:
        v = table[n]
        print(n, v["R180C"], v["R180M"], v["R180V"], v["total"])

    for col, (desc, oeis, parity) in COLUMNS.items():
        name = f"{col} ({oeis})" if oeis else f"{col} (no OEIS number known)"
        print()
        print(f"# {name}: {desc}, n = {ns[0]}..{ns[-1]}")
        print(", ".join(str(table[n][col]) for n in ns))
        if parity is not None:
            sub = [n for n in ns if n % 2 == parity]
            for n in ns:
                if n % 2 != parity and table[n][col] != 0:
                    sys.exit(f"{col} is nonzero for n={n}?")
            print(f"# {col}, only {'odd' if parity else 'even'} n = {sub[0]}, {sub[1]}, ..., {sub[-1]}")
            print(", ".join(str(table[n][col]) for n in sub))
        if prefix:
            bname = f"{prefix}-{oeis or col}.txt"
            with open(bname, "w") as f:
                for n in ns:
                    f.write(f"{n} {table[n][col]}\n")
            print(f"# b-file written: {bname}")


if __name__ == "__main__":
    main()
