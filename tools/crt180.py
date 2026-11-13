"""Combine diag180 results modulo several moduli (file results180.py) by CRT.

usage: python3 tools/crt180.py [results180.py] [bits] [-b PREFIX] [-c DIR]

Prints
  1. a table n, R180C, R180M, R180V, total;
  2. for each OEIS sequence (and for the total) the bare comma-separated
     list of numbers, as wanted for the OEIS data field;
  3. with -b PREFIX, b-files PREFIXb390653.txt, ... in OEIS format
     (lines "n a(n)");
  4. with -c DIR, a comparison with the OEIS b-files b390653.txt, ... in the
     directory DIR (the repository has them in b-files/).

A value is accepted only if its CRT reconstruction is below 2**bits
(default 72); a wrong run would give a value of the size of the product of
the moduli.  A result file with modulus 2**128 (128-bit build) needs no
second modulus.

The columns count FIXED polyominoes that are invariant under the rotation by
180 degrees (whatever other symmetries they have), by the type of center.
They are in the OEIS (John Mason, 2025, from Shirakawa's table):
    A390653(n) = R180C(n)        center in a cell
    A390654(n) = R180M(2n) / 2   center on the midpoint of a horizontal edge
                                 (R180M counts horizontal and vertical edges)
    A390655(n) = R180V(2n)       center at a vertex
"""
import sys
from math import prod

# name, OEIS number, description, index n -> (polyomino size, column, divisor)
SEQUENCES = [
    ("R180C", "A390653", "n-celled fixed polyominoes with at least 180 degree rotational symmetry about the center of a square",
     lambda n: (n, "R180C", 1)),
    ("R180M/2", "A390654", "2n-celled fixed polyominoes with at least 180 degree rotational symmetry about the midpoint of a horizontal edge",
     lambda n: (2 * n, "R180M", 2)),
    ("R180V", "A390655", "2n-celled fixed polyominoes with at least 180 degree rotational symmetry about a vertex",
     lambda n: (2 * n, "R180V", 1)),
    ("total", None, "n-celled fixed polyominoes with at least 180 degree rotational symmetry",
     lambda n: (n, "total", 1)),
]


def read_bfile(fn):
    d = {}
    for line in open(fn):
        if line.strip() and not line.startswith("#"):
            n, v = line.split()[:2]
            d[int(n)] = int(v)
    return d


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
    prefix = checkdir = None
    for opt in ("-b", "-c"):
        if opt in args:
            i = args.index(opt)
            if opt == "-b":
                prefix = args[i + 1]
            else:
                checkdir = args[i + 1]
            del args[i:i + 2]
    fn = args[0] if len(args) > 0 else "results180.py"
    bits = int(args[1]) if len(args) > 1 else 72
    ns, mods, table = combine(fn, bits)

    print("# moduli:", ", ".join(str(m) for m in mods))
    print("# n R180C R180M R180V total")
    for n in ns:
        v = table[n]
        print(n, v["R180C"], v["R180M"], v["R180V"], v["total"])

    for n in ns:
        if n % 2 and (table[n]["R180M"] or table[n]["R180V"]):
            sys.exit(f"R180M or R180V is nonzero for odd n={n}?")
        if table[n]["R180M"] % 2:
            sys.exit(f"R180M({n}) is odd?")
    for name, oeis, desc, idx in SEQUENCES:
        seq = {}
        k = 1
        while idx(k)[0] <= ns[-1]:
            n, col, div = idx(k)
            if n in table:
                seq[k] = table[n][col] // div
            k += 1
        ks = sorted(seq)
        print()
        print(f"# {oeis or 'no OEIS number'}: {name}, {desc}; n = {ks[0]}..{ks[-1]}")
        print(", ".join(str(seq[k]) for k in ks))
        if prefix is not None:
            bname = f"{prefix}b{oeis[1:]}.txt" if oeis else f"{prefix}b-{name}.txt"
            with open(bname, "w") as f:
                for k in ks:
                    f.write(f"{k} {seq[k]}\n")
            print(f"# b-file written: {bname}")
        if checkdir is not None and oeis:
            old = read_bfile(f"{checkdir}/b{oeis[1:]}.txt")
            common = [k for k in ks if k in old]
            bad = [k for k in common if old[k] != seq[k]]
            new = [k for k in ks if k not in old]
            print(f"# compared with {checkdir}/b{oeis[1:]}.txt: {len(common)} common terms, "
                  + (f"MISMATCH at n = {bad}" if bad else "all equal")
                  + (f"; new terms n = {new[0]}..{new[-1]}" if new else ""))


if __name__ == "__main__":
    main()
