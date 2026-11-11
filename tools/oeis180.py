"""Derive OEIS sequences for FREE polyominoes with exactly C2 symmetry from the
fixed counts of diag180, and check them against the OEIS b-files.

usage: python3 tools/oeis180.py NEWDIR/ [BFILEDIR]

NEWDIR/ is the prefix given to crt180.py -b (it contains b390653.txt,
b390654.txt, b390655.txt); BFILEDIR (default b-files) has the OEIS b-files.

Formulas (John Mason, OEIS):
  A390653(n) = 4*A351615(n) + 2*A351190(n) + 2*A351159(n) + 2*A351142(n) + A351127(n)
  A390654(n) = 2*A234008(n) + A346799(n)
  A390655(n) = 4*A351616(n) + 2*A351191(n/2) + 2*A351160(n) + 2*A234007(n/2) + A346800(n/2)
  A006747(n) = A351615(n) + A234008(n/2) + A351616(n/2)   (last two for even n)
The divisions by 4 and 2 must be exact, and the old terms must be
reproduced: both are checked.  New terms are printed.
"""
import sys


def bf(fn):
    d = {}
    for line in open(fn):
        p = line.split()
        if len(p) >= 2 and not line.startswith("#"):
            d[int(p[0])] = int(p[1])
    return d


new = sys.argv[1]
old = sys.argv[2] if len(sys.argv) > 2 else "b-files"
C, M, V = (bf(f"{new}b{a}.txt") for a in ("390653", "390654", "390655"))
B = {a: bf(f"{old}/b{a}.txt") for a in
     "351615 351190 351159 351142 351127 351616 351191 351160 234007 346800 346799 006747".split()}


def half(a, n):
    return B[a].get(n // 2) if n % 2 == 0 else 0


def derive(name, src, parts_of, div):
    out, checked = {}, 0
    for n in sorted(src):
        parts = parts_of(n)
        if None in parts:
            break
        x = src[n] - sum(parts)
        if x % div:
            sys.exit(f"{name}({n}): not divisible by {div}")
        out[n] = x // div
        if name in B and n in B[name]:
            if B[name][n] != out[n]:
                sys.exit(f"A{name}({n}) differs from the OEIS b-file")
            checked += 1
    newterms = {n: v for n, v in out.items() if name not in B or n not in B[name]}
    if name not in B:
        print(f"A{name}: {len(out)} terms derived (no b-file to compare with):")
    else:
        print(f"A{name}: {len(out)} terms derived, {checked} equal to the b-file; "
              + (f"terms not in the b-file: n = {min(newterms)}..{max(newterms)}" if newterms else "no new terms"))
    for n, v in newterms.items():
        print(f"  {n} {v}")
    return out


a615 = derive("351615", C, lambda n: [2 * B["351190"].get(n) if n in B["351190"] else None,
                                      2 * B["351159"].get(n) if n in B["351159"] else None,
                                      2 * B["351142"].get(n) if n in B["351142"] else None,
                                      B["351127"].get(n)], 4)
a008 = derive("234008", M, lambda n: [B["346799"].get(n)], 2)


def vparts(n):
    p = [half("351191", n), B["351160"].get(n), half("234007", n), half("346800", n)]
    if None in p:
        return [None]
    return [2 * p[0], 2 * p[1], 2 * p[2], p[3]]


a616 = derive("351616", V, vparts, 4)
tot = {}
for n in sorted(a615):
    if n % 2 == 0:
        if n // 2 not in a008 or n // 2 not in a616:
            break
        tot[n] = a615[n] + a008[n // 2] + a616[n // 2]
    else:
        tot[n] = a615[n]
bad = [n for n in tot if n in B["006747"] and B["006747"][n] != tot[n]]
if bad:
    sys.exit(f"A006747 differs from the b-file at n = {bad}")
newt = {n: v for n, v in tot.items() if n not in B["006747"]}
print(f"A006747: {len(tot)} terms derived, {len(tot) - len(newt)} equal to the b-file; new terms:")
for n, v in newt.items():
    print(f"  {n} {v}")
