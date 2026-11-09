/*
 * redelmeier.c -- baseline polyomino counter (Redelmeier's algorithm).
 *
 * Counts fixed polyominoes (OEIS A001168) of every size up to N.
 * With -s it also counts, for each size, how many fixed polyominoes are
 * invariant under each element of the symmetry group D4 of the square, and
 * derives the number of free polyominoes (OEIS A000105) by Burnside's lemma:
 *
 *   free = (fix(id) + fix(r90) + fix(r180) + fix(r270)
 *           + fix(mirror_v) + fix(mirror_h) + fix(diag) + fix(antidiag)) / 8
 *
 * where fix(r90) = fix(r270).  This symmetry test runs on every polyomino
 * the enumeration produces, so -s is a slow reference implementation; it is
 * meant as ground truth for faster methods, not as the production counter.
 *
 * Usage: redelmeier [-s] N
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#define MAXN 40

typedef unsigned long long u64;

static int N;            /* maximum size */
static int W;            /* row width of the board array */
static int do_sym;

static unsigned char *mark;   /* 1 = cell is (or was) in the untried set, or forbidden */
static int poly[MAXN];        /* cells of the current polyomino */
static int *untried[MAXN + 1];/* untried set per level */

static u64 count[MAXN + 1];

/* symmetry statistics, per size */
enum { S_R90, S_R180, S_MV, S_MH, S_D, S_AD, NSYM };
static const char *symname[NSYM] = { "r90", "r180", "mirror_v", "mirror_h", "diag", "antidiag" };
static u64 symcount[MAXN + 1][NSYM];
/* fixed polyominoes by their exact set of symmetries (bitmask over NSYM) */
static u64 groupcount[MAXN + 1][1 << NSYM];

static unsigned char *grid;   /* N x N scratch bitmap for symmetry tests */

static void check_symmetry(int k)
{
    int xs[MAXN], ys[MAXN];
    int minx = 1 << 30, maxx = -1, miny = 1 << 30, maxy = -1;
    for (int i = 0; i < k; i++) {
        int x = poly[i] % W, y = poly[i] / W;
        xs[i] = x; ys[i] = y;
        if (x < minx) minx = x;
        if (x > maxx) maxx = x;
        if (y < miny) miny = y;
        if (y > maxy) maxy = y;
    }
    int w = maxx - minx + 1, h = maxy - miny + 1;
    for (int i = 0; i < k; i++) {
        xs[i] -= minx; ys[i] -= miny;
        grid[ys[i] * N + xs[i]] = 1;
    }
#define CELL(x, y) grid[(y) * N + (x)]
    int ok[NSYM];
    for (int s = 0; s < NSYM; s++) ok[s] = 1;
    if (w != h) ok[S_R90] = ok[S_D] = ok[S_AD] = 0;
    for (int i = 0; i < k; i++) {
        int x = xs[i], y = ys[i];
        if (ok[S_R180] && !CELL(w - 1 - x, h - 1 - y)) ok[S_R180] = 0;
        if (ok[S_MV]   && !CELL(w - 1 - x, y))         ok[S_MV] = 0;
        if (ok[S_MH]   && !CELL(x, h - 1 - y))         ok[S_MH] = 0;
        if (ok[S_R90]  && !CELL(y, w - 1 - x))         ok[S_R90] = 0;
        if (ok[S_D]    && !CELL(y, x))                 ok[S_D] = 0;
        if (ok[S_AD]   && !CELL(w - 1 - y, w - 1 - x)) ok[S_AD] = 0;
    }
#undef CELL
    int mask = 0;
    for (int s = 0; s < NSYM; s++)
        if (ok[s]) { symcount[k][s]++; mask |= 1 << s; }
    groupcount[k][mask]++;
    for (int i = 0; i < k; i++) grid[ys[i] * N + xs[i]] = 0;
}

static const int *nbr_off;    /* the 4 neighbour offsets */

/* depth = number of cells already in the polyomino; u = untried set, nu = its size */
static void rec(int depth, int *u, int nu)
{
    int *nu_set = untried[depth + 1];
    while (nu > 0) {
        int c = u[--nu];              /* untried set used as a stack */
        poly[depth] = c;
        int k = depth + 1;
        count[k]++;
        if (do_sym) check_symmetry(k);
        if (k < N) {
            memcpy(nu_set, u, nu * sizeof(int));
            int m = nu;
            for (int j = 0; j < 4; j++) {
                int d = c + nbr_off[j];
                if (!mark[d]) { mark[d] = 1; nu_set[m++] = d; }
            }
            rec(k, nu_set, m);
            for (int j = nu; j < m; j++) mark[nu_set[j]] = 0;
        }
    }
}

int main(int argc, char **argv)
{
    int a = 1;
    if (a < argc && strcmp(argv[a], "-s") == 0) { do_sym = 1; a++; }
    if (a >= argc) { fprintf(stderr, "usage: %s [-s] N\n", argv[0]); return 1; }
    N = atoi(argv[a]);
    if (N < 1 || N > MAXN) { fprintf(stderr, "N must be in 1..%d\n", MAXN); return 1; }

    /* Board: x in [-N, N] stored as x+N (width 2N+1), y in [-1, N] stored as y+1.
       Cells with y < 0, or y == 0 and x < 0, are forbidden: the origin (0,0)
       is the lowest-then-leftmost cell of every polyomino, which makes each
       fixed polyomino appear exactly once. */
    W = 2 * N + 1;
    int H = N + 2;
    mark = calloc((size_t)W * H, 1);
    grid = calloc((size_t)N * N, 1);
    for (int i = 0; i <= N; i++) untried[i] = malloc(sizeof(int) * (4 * N + 4));
    static int off[4];
    off[0] = 1; off[1] = -1; off[2] = W; off[3] = -W;
    nbr_off = off;
    for (int x = 0; x < W; x++) mark[x] = 1;          /* row y = -1 */
    for (int x = 0; x < N; x++) mark[W + x] = 1;      /* y = 0, x < 0 */
    int origin = W + N;
    mark[origin] = 1;
    untried[0][0] = origin;

    clock_t t0 = clock();
    rec(0, untried[0], 1);
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;

    if (!do_sym) {
        printf("# n fixed\n");
        for (int k = 1; k <= N; k++) printf("%d %llu\n", k, count[k]);
    } else {
        printf("# n fixed");
        for (int s = 0; s < NSYM; s++) printf(" %s", symname[s]);
        printf(" free\n");
        int bad = 0;
        for (int k = 1; k <= N; k++) {
            u64 *f = symcount[k];
            u64 sum = count[k] + 2 * f[S_R90] + f[S_R180] + f[S_MV] + f[S_MH] + f[S_D] + f[S_AD];
            printf("%d %llu", k, count[k]);
            for (int s = 0; s < NSYM; s++) printf(" %llu", f[s]);
            printf(" %llu%s\n", sum / 8, sum % 8 ? " (NOT DIVISIBLE BY 8!)" : "");
            if (sum % 8) bad = 1;
        }
        /* Free polyominoes by symmetry type.  A free polyomino whose symmetry
           group has order g has 8/g fixed versions; for the types with a
           single mirror axis (or a pair of axes) these versions are split
           between conjugate masks, so masks are merged per type first. */
        static const struct { const char *name; int m1, m2, order; } types[] = {
            { "none",       0, -1, 1 },
            { "r180",       1 << S_R180, -1, 2 },
            { "r90",        (1 << S_R90) | (1 << S_R180), -1, 4 },
            { "mirror_axis", 1 << S_MV, 1 << S_MH, 2 },
            { "mirror_diag", 1 << S_D, 1 << S_AD, 2 },
            { "two_axes",   (1 << S_R180) | (1 << S_MV) | (1 << S_MH), -1, 4 },
            { "two_diags",  (1 << S_R180) | (1 << S_D) | (1 << S_AD), -1, 4 },
            { "all",        (1 << NSYM) - 1, -1, 8 },
        };
        int ntypes = sizeof types / sizeof types[0];
        printf("# free polyominoes by symmetry type\n# n");
        for (int t = 0; t < ntypes; t++) printf(" %s", types[t].name);
        printf("\n");
        for (int k = 1; k <= N; k++) {
            u64 total = 0, seen = 0;
            printf("%d", k);
            for (int t = 0; t < ntypes; t++) {
                u64 c = groupcount[k][types[t].m1];
                if (types[t].m2 >= 0) c += groupcount[k][types[t].m2];
                seen += c;
                if (c * types[t].order % 8) bad = 1;
                c = c * types[t].order / 8;
                total += c;
                printf(" %llu", c);
            }
            printf("\n");
            if (seen != count[k]) { fprintf(stderr, "unexpected symmetry mask at n=%d\n", k); bad = 1; }
        }
        if (bad) { fprintf(stderr, "consistency check FAILED\n"); return 2; }
    }
    fprintf(stderr, "time %.2f s\n", secs);
    return 0;
}
