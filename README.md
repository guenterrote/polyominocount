# polyominocount
counting polyominos, with the goal to figure out numbers of free polyominos

## Baseline: Redelmeier's algorithm (`src/redelmeier.c`)

```
make            # builds build/redelmeier
make test       # compares with a brute-force Python reference (tests/brute.py)
build/redelmeier 18      # fixed polyominoes (A001168) for n = 1..18
build/redelmeier -s 14   # + symmetry counts, free polyominoes via Burnside
```

The untried set is used as a stack. With `-s`, every polyomino produced by
the enumeration is tested against the 7 non-trivial elements of D4. This
yields, for each n, the number of fixed polyominoes invariant under r90,
r180, the two axis mirrors and the two diagonal mirrors, then

    free = (fixed + 2*r90 + r180 + mirror_v + mirror_h + diag + antidiag) / 8

and a breakdown of free polyominoes by symmetry type. Divisibility by 8 and
the type breakdown are checked as a consistency test. The `-s` mode is a
slow reference (ground truth for faster methods), not a production counter:
symmetric polyominoes are an exponentially small fraction and should be
counted directly.

### Timings (single thread, cloud container, gcc -O3)

| n  | fixed only | with `-s` |
|----|-----------:|----------:|
| 14 | 0.06 s     | 1.4 s     |
| 16 | 0.85 s     | 19 s      |
| 18 | 11.7 s     |           |
| 19 | 46.5 s     |           |
| 20 | 176 s      |           |

Growth is about 3.8x per extra cell, i.e. roughly 1.6e8 polyominoes per second.

Verified: fixed counts match A001168 for n <= 20; free counts match A000105
for n <= 12 and the brute-force reference for n <= 11.

## State of the art (checked Nov 2026)

- Fixed polyominoes A001168: Jensen (transfer matrix) up to n = 56,
  Barequet and Ben-Shachar extended to n = 70.
- Free polyominoes A000105: known up to n = 59 (T. Shirakawa,
  arXiv:2510.22446, Oct 2025), via the symmetry-class formula of John Mason:
  a(n) = (4*A006746 + 4*A006748 + 4*A006747 + 6*A056877 + 6*A056878
  + 6*A144553 + 7*A142886 + A001168)/8.
- So for 60 <= n <= 70 the fixed counts exist and only the symmetric classes
  are missing.

## diag180: polyominoes with 180-degree symmetry (CWEB)

`src/diag180.w` counts the fixed polyominoes that are invariant under a
rotation by 180 degrees, with a transfer matrix along diagonal rows. It is
written in CWEB with the LaTeX class cweb-hy (files in `cweb/`), like
[NumPSLA](https://github.com/guenterrote/NumPSLA). The typeset program with
explanations, figures and results is `doc/diag180.pdf`.

- `make build/diag180` compiles `src/diag180.c`, which is the output of
  ctangle and is kept in the repository. `make CNT128=1 build/diag180` uses
  128-bit counters.
- `make tangle` regenerates `src/diag180.c` from `src/diag180.w`.
- `make doc` regenerates `doc/diag180.pdf` and `doc/diag180p.pdf` (cweave,
  pdflatex, bibtex with the style plainurl).

### Two versions from one source

`src/production.ch` is a CWEB change file. Together with `src/diag180.w` it
gives the production version `diag180p` (`src/diag180p.c`,
`doc/diag180p.pdf`, `make build/diag180p`), which lacks the experimental
options `-K`, `-A`, `-N` and `-H`. There is only one text: an improvement in
`diag180.w` appears in both versions after `make tangle` and `make doc`.
If a passage that the change file replaces is edited in `diag180.w`, ctangle
and cweave stop with "Change file entry did not match"; then copy the new
lines into the corresponding `@x` part of `production.ch`. `make test` checks
that both versions give the same counts.
- `tools/crt180.py` combines the result files of runs with different moduli.
