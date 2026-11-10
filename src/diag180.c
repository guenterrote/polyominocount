/*
 * diag180: count polyominoes with 180-degree rotational symmetry, by an
 * outward transfer matrix along DIAGONAL rows.  C version of
 * proto/rot180_diag.py (same states, same pruning, same results).
 *
 * Coordinates u = x + y, v = x - y.  Neighbours of (u, v) are (u+-1, v+-1).
 * rho(x, y) = (c - x, d - y) becomes (u, v) -> (s - u, t - v), s = c + d,
 * t = c - d.  The strip is v in [t - vhi, vhi], L = 2 vhi - t + 1 positions,
 * P touches v = vhi (and by symmetry v = t - vhi).  We count the half
 * H = cells with u > s/2 (plus the middle row if s is even) on the
 * quotient P / rho: every frontier class has a partner rho(class), which is
 * itself ('s'), another frontier class, or a ghost ('g') without frontier
 * cells.
 *
 * fix(r180)(n) = sum over (c, d, vhi) of 2 #(vspan < uspan) + #(vspan = uspan).
 *
 * Counters are uint64 modulo 2^64, or modulo M with -m M (for CRT).
 *
 * usage: diag180 [-m M] [-R] N              all cases, prints n and fix(r180)(n)
 *        diag180 [-m M] [-R] N c d vhi      one case, prints n, #lt, #eq
 *   -R  disable the reflection merge (v -> -v) for t = 0
 *   -N  no one-step lookahead in the bound (default: lookahead at the kink)
 *   -A  one-step lookahead after every cell
 *   -K  use the older, weaker bound (MST only, no self/kink terms)
 *   -C  check that every packed key unpacks to the same state
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

typedef uint64_t cnt_t;

static int NMAX;
static uint64_t MOD = 0;          /* 0: arithmetic mod 2^64 */
static int REFLECT = 1;
static int OLDBOUND = 0;
static int LOOKAHEAD = 1;       /* 1: lookahead at the kink, 2 (-A): everywhere, 0 (-N): none */        /* -K: Kruskal bound without self/kink terms */

static inline cnt_t addm(cnt_t a, cnt_t b)
{
    if (!MOD)
        return a + b;
    cnt_t s = a + b;
    if (s < a || s >= MOD)
        s -= MOD;
    return s;
}

/* ------------------------------------------------------------------ */
/* case geometry */
static int L, vlo, T, S, E;       /* E = s mod 2 */
static int MAXC;                  /* max classes in a state */
static int KL;                    /* key length in bytes */
static int NW;                    /* counter vector length = NMAX + 1 */

#define MAXL 160
#define PSELF 0xFF
#define PGHOST 0xFE

/* ------------------------------------------------------------------ */
/* State storage.  A canonical state is first built as a byte key (L labels,
   MAXC partner codes, touched flag).  In the table it is bit-packed:
     per position 3 bits: 0 empty, 1 continues the run to its left,
       run starts: 2 class with a single run, 3 first run of a class,
       4 middle run, 5 last run (classes are non-crossing, so a stack
       recovers the labels);
     per class 1 bit: partner present (or self) / ghost.  Classes with a
       present partner or self-partnered cross the disk between the
       frontier and its rho-image, and rho reverses their order, so they
       form one nested chain A1 < ... < Ak < (S) < rho Ak < ... < rho A1:
       the i-th crossing class from the left is paired with the i-th from
       the right, and the middle one (odd count) is self;
     1 bit touched.
   Counts per state are a chain of (size, count) records. */
static int KW;                    /* packed key length in 64-bit words */

typedef struct {
    cnt_t c;
    uint32_t next;                /* index + 1, 0 = end */
    uint8_t wt;
} rec_t;

typedef struct {
    size_t n, cap, hmask;
    uint64_t *keys;
    uint32_t *head;               /* first record + 1 */
    rec_t *rec;
    size_t nrec, rcap;
    uint32_t *hash;               /* index + 1, 0 = empty */
} table_t;

static void *xrealloc(void *p, size_t s)
{
    p = realloc(p, s);
    if (!p) {
        fprintf(stderr, "out of memory\n");
        exit(1);
    }
    return p;
}

static void tab_init(table_t *t, size_t cap)
{
    t->n = 0;
    t->cap = cap;
    size_t h = 16;
    while (h < 2 * cap)
        h <<= 1;
    t->hmask = h - 1;
    t->keys = xrealloc(0, cap * KW * 8);
    t->head = xrealloc(0, cap * 4);
    t->hash = calloc(h, sizeof(uint32_t));
    t->rcap = 2 * cap;
    t->nrec = 0;
    t->rec = xrealloc(0, t->rcap * sizeof(rec_t));
}

static void tab_free(table_t *t)
{
    free(t->keys);
    free(t->head);
    free(t->rec);
    free(t->hash);
}

static size_t tab_bytes(const table_t *t)
{
    return t->n * (KW * 8 + 4) + t->nrec * sizeof(rec_t) + (t->hmask + 1) * 4;
}

static inline uint64_t hkey(const uint64_t *k)
{
    uint64_t h = 0x9E3779B97F4A7C15ull;
    for (int i = 0; i < KW; i++) {
        h ^= k[i];
        h *= 0xBF58476D1CE4E5B9ull;
        h ^= h >> 31;
    }
    return h;
}

static void tab_grow(table_t *t)
{
    size_t cap = t->cap * 2;
    t->keys = xrealloc(t->keys, cap * KW * 8);
    t->head = xrealloc(t->head, cap * 4);
    size_t h = (t->hmask + 1) * 2;
    free(t->hash);
    t->hash = calloc(h, sizeof(uint32_t));
    if (!t->hash) {
        fprintf(stderr, "out of memory\n");
        exit(1);
    }
    t->hmask = h - 1;
    t->cap = cap;
    for (size_t i = 0; i < t->n; i++) {
        size_t p = hkey(t->keys + i * KW) & t->hmask;
        while (t->hash[p])
            p = (p + 1) & t->hmask;
        t->hash[p] = i + 1;
    }
}

static void tab_clear(table_t *t)
{
    t->n = 0;
    t->nrec = 0;
    memset(t->hash, 0, (t->hmask + 1) * sizeof(uint32_t));
}

static inline void putbits(uint64_t *w, int *pos, uint64_t v, int nb)
{
    int i = *pos >> 6, o = *pos & 63;
    w[i] |= v << o;
    if (o + nb > 64)
        w[i + 1] |= v >> (64 - o);
    *pos += nb;
}

static inline int getbits(const uint64_t *w, int *pos, int nb)
{
    int i = *pos >> 6, o = *pos & 63;
    uint64_t v = w[i] >> o;
    if (o + nb > 64)
        v |= w[i + 1] << (64 - o);
    *pos += nb;
    return v & ((1u << nb) - 1);
}

static void pack(const uint8_t *key, uint64_t *w)
{
    int first[MAXL + 2], lastrun[MAXL + 2];
    int pos = 0, m = 0;
    memset(w, 0, KW * 8);
    for (int i = 0; i < L; i++)
        if (key[i] && (i == 0 || !key[i - 1])) {
            int l = key[i];
            if (l > m) {
                m = l;
                first[l] = i;
            }
            lastrun[l] = i;
        }
    for (int i = 0; i < L; i++) {
        int l = key[i], code;
        if (!l)
            code = 0;
        else if (i > 0 && key[i - 1])
            code = 1;
        else if (first[l] == i)
            code = lastrun[l] == i ? 2 : 3;
        else
            code = lastrun[l] == i ? 5 : 4;
        putbits(w, &pos, code, 3);
    }
    for (int j = 1; j <= m; j++)
        putbits(w, &pos, key[L + j - 1] != PGHOST, 1);
    putbits(w, &pos, key[KL - 1], 1);
}

static void unpack(const uint64_t *w, uint8_t *key)
{
    int stack[MAXL + 2], sp = 0, m = 0, pos = 0;
    memset(key, 0, KL);
    for (int i = 0; i < L; i++) {
        int code = getbits(w, &pos, 3);
        switch (code) {
        case 0: key[i] = 0; break;
        case 1: key[i] = key[i - 1]; break;
        case 2: key[i] = ++m; break;
        case 3: key[i] = ++m; stack[sp++] = m; break;
        case 4: key[i] = stack[sp - 1]; break;
        case 5: key[i] = stack[--sp]; break;
        }
    }
    sp = 0;                       /* crossing classes, left to right */
    for (int j = 1; j <= m; j++)
        if (getbits(w, &pos, 1))
            stack[sp++] = j;
        else
            key[L + j - 1] = PGHOST;
    for (int a = 0, b = sp - 1; a <= b; a++, b--) {
        key[L + stack[a] - 1] = a == b ? PSELF : stack[b];
        key[L + stack[b] - 1] = a == b ? PSELF : stack[a];
    }
    key[KL - 1] = getbits(w, &pos, 1);
}

static int CHECK = 0;

/* index of the state with this byte key (created if new) */
static size_t tab_get(table_t *t, const uint8_t *key)
{
    uint64_t w[12];
    pack(key, w);
    if (CHECK) {
        uint8_t k2[2 * MAXL + 4];
        unpack(w, k2);
        if (memcmp(k2, key, KL)) {
            fprintf(stderr, "pack/unpack mismatch\n");
            exit(2);
        }
    }
    size_t p = hkey(w) & t->hmask;
    while (t->hash[p]) {
        size_t i = t->hash[p] - 1;
        if (!memcmp(t->keys + i * KW, w, KW * 8))
            return i;
        p = (p + 1) & t->hmask;
    }
    if (t->n == t->cap) {
        tab_grow(t);
        p = hkey(w) & t->hmask;
        while (t->hash[p])
            p = (p + 1) & t->hmask;
    }
    size_t i = t->n++;
    t->hash[p] = i + 1;
    memcpy(t->keys + i * KW, w, KW * 8);
    t->head[i] = 0;
    return i;
}

static void tab_add(table_t *t, size_t si, int wt, cnt_t c)
{
    for (uint32_t r = t->head[si]; r; r = t->rec[r - 1].next)
        if (t->rec[r - 1].wt == wt) {
            t->rec[r - 1].c = addm(t->rec[r - 1].c, c);
            return;
        }
    if (t->nrec == t->rcap) {
        t->rcap *= 2;
        t->rec = xrealloc(t->rec, t->rcap * sizeof(rec_t));
    }
    rec_t *r = &t->rec[t->nrec++];
    r->c = c;
    r->wt = wt;
    r->next = t->head[si];
    t->head[si] = t->nrec;
}

/* dense counter vector of state si (gf must be zero on entry) */
static void tab_load(const table_t *t, size_t si, cnt_t *gf)
{
    for (uint32_t r = t->head[si]; r; r = t->rec[r - 1].next)
        gf[t->rec[r - 1].wt] = t->rec[r - 1].c;
}

static void tab_unload(const table_t *t, size_t si, cnt_t *gf)
{
    for (uint32_t r = t->head[si]; r; r = t->rec[r - 1].next)
        gf[t->rec[r - 1].wt] = 0;
}

/* ------------------------------------------------------------------ */
/* working state: frontier labels, partners, union-find over classes and
   ghosts.  Nodes 1..MAXC+1 are classes, GH(i) = ghost of class i. */
#define MAXN 2 * (MAXL + 4)
static int GOFF;                  /* ghost(i) = GOFF + i */

typedef struct {
    int f[MAXL];
    int partner[MAXL];            /* for classes: -1 self, -2 ghost, else class */
    int parent[MAXN];
    int next;                     /* next free class label */
} work_t;

#define GH(i) (GOFF + (i))
#define ISGHOST(a) ((a) > GOFF)

static inline int find(work_t *w, int a)
{
    while (w->parent[a] != a) {
        w->parent[a] = w->parent[w->parent[a]];
        a = w->parent[a];
    }
    return a;
}

static inline int rho(work_t *w, int a)
{
    a = find(w, a);
    if (ISGHOST(a))
        return find(w, a - GOFF);
    int p = w->partner[a];
    if (p == -1)
        return a;
    if (p == -2)
        return find(w, GH(a));
    return find(w, p);
}

static inline void link_(work_t *w, int a, int b)
{
    a = find(w, a);
    b = find(w, b);
    if (a == b)
        return;
    if (ISGHOST(a)) {
        int x = a;
        a = b;
        b = x;
    }
    w->parent[b] = a;
}

static inline void unite(work_t *w, int a, int b)
{
    int ra = rho(w, a), rb = rho(w, b);
    link_(w, a, b);
    link_(w, ra, rb);
}

static void decode(work_t *w, const uint8_t *key)
{
    int m = 0;
    for (int i = 0; i < L; i++) {
        w->f[i] = key[i];
        if (key[i] > m)
            m = key[i];
    }
    for (int a = 0; a < MAXN; a++)
        w->parent[a] = a;
    for (int j = 1; j <= MAXC; j++) {
        int q = key[L + j - 1];
        w->partner[j] = q == PSELF ? -1 : q == PGHOST ? -2 : q;
        if (q && j > m)
            m = j;
    }
    w->next = m + 1;
}

/* canonical key; returns number of classes */
static int encode(work_t *w, uint8_t *key, int tch)
{
    int map[MAXN];
    int rep[MAXL];
    int m = 0;
    memset(map, 0, sizeof(map));
    memset(key, 0, KL);
    for (int i = 0; i < L; i++) {
        if (!w->f[i])
            continue;
        int r = find(w, w->f[i]);
        if (!map[r]) {
            map[r] = ++m;
            rep[m] = r;
        }
        key[i] = map[r];
    }
    for (int j = 1; j <= m; j++) {
        int q = rho(w, rep[j]);
        key[L + j - 1] = q == rep[j] ? PSELF : map[q] ? map[q] : PGHOST;
    }
    key[KL - 1] = tch;
    return m;
}

enum { OK, DEAD, COMPLETE };

static int retire(work_t *w, int i)
{
    int lab = w->f[i];
    if (!lab)
        return OK;
    w->f[i] = 0;
    int r = find(w, lab);
    int any = 0;
    for (int j = 0; j < L; j++)
        if (w->f[j]) {
            any = 1;
            if (find(w, w->f[j]) == r)
                return OK;
        }
    int p = rho(w, r);
    if (p == r)
        return any ? DEAD : COMPLETE;
    if (ISGHOST(p))
        return DEAD;
    w->partner[find(w, p)] = -2;
    return OK;
}

/* ------------------------------------------------------------------ */
/* bounds on a canonical key */
static int completion_possible(const uint8_t *key)
{
    int m = 0;
    for (int j = 0; j < MAXC; j++) {
        int q = key[L + j];
        if (!q)
            break;
        m++;
        if (q != PGHOST)
            return 1;
    }
    return m == 0;
}

/* exact number of further H cells needed to join the super-classes and
   (if not yet touched) to reach an edge of the strip */
static int exact_bound(const uint8_t *key, int tch)
{
    int sup[MAXL + 1];
    int occ[MAXL], no = 0;
    for (int j = 1; j <= MAXC; j++) {
        int q = key[L + j - 1];
        if (!q)
            break;
        sup[j] = (q != PSELF && q != PGHOST && q < j) ? q : j;
    }
    for (int i = 0; i < L; i++)
        if (key[i])
            occ[no++] = i;
    if (!no)
        return 0;
    /* Kruskal over gaps between consecutive occupied positions */
    int par[MAXL + 1];
    int ncls = 0;
    int seen[MAXL + 1];
    memset(seen, 0, sizeof(seen));
    for (int a = 0; a < no; a++) {
        int s = sup[key[occ[a]]];
        if (!seen[s]) {
            seen[s] = 1;
            ncls++;
        }
        par[s] = s;
    }
    int cost = 0;
    if (ncls > 1) {
        /* counting sort of gaps by weight */
        int ne = no - 1;
        int wgt[MAXL], ord[MAXL], bucket[MAXL + 1];
        memset(bucket, 0, sizeof(bucket));
        for (int a = 0; a < ne; a++) {
            wgt[a] = occ[a + 1] - occ[a] - 1;
            bucket[wgt[a]]++;
        }
        for (int g = 1; g <= L; g++)
            bucket[g] += bucket[g - 1];
        for (int a = ne - 1; a >= 0; a--)
            ord[--bucket[wgt[a]]] = a;
        for (int e = 0; e < ne && ncls > 1; e++) {
            int a = ord[e];
            int x = sup[key[occ[a]]], y = sup[key[occ[a + 1]]];
            while (par[x] != x)
                x = par[x];
            while (par[y] != y)
                y = par[y];
            if (x != y) {
                par[x] = y;
                cost += wgt[a];
                ncls--;
            }
        }
    }
    if (!tch) {
        int d1 = occ[0], d2 = L - 1 - occ[no - 1];
        cost += d1 < d2 ? d1 : d2;
    }
    return cost;
}

/* Cells of H still needed (Rote): on the diagonal lattice, connections
   between frontier runs use disjoint v columns, one cell per column.  The
   graph of gaps between consecutive runs (classes identified with their
   partners) is a cactus, because classes are non-crossing:
     - the pocket between two consecutive runs of a class X, containing
       classes C1..Cr, is the cycle X - C1 - ... - Cr - X: pay the sum of
       its gaps minus the largest;
     - at the top level, the crossing classes nest as
       A1 < ... < Ak < (S) < rho Ak < ... < rho A1; outer gaps (left of A1,
       right of rho A1) are paths and are paid in full; ring i (A_i..A_i+1
       and rho A_i+1..rho A_i) and the middle (A_k..S..rho A_k, or
       A_k..rho A_k) are cycles: sum minus max;
     - without a self class some X must be joined to rho X, which closes
       one of the top-level cycles completely: add the smallest of their
       maxima.
   kink >= 0: the gap ending at position kink costs one more cell (the row
   u-1 cell there is reachable from the left only over its own column).
   Plus, if the strip edge is not yet touched, the distance to it. */
static int cactus_bound(const uint8_t *key, int tch, int kink)
{
    int rs[MAXL], re[MAXL], rl[MAXL], R = 0;
    int lastrun[MAXL + 2];
    for (int i = 0; i < L; i++)
        if (key[i]) {
            if (i == 0 || !key[i - 1]) {
                rs[R] = i;
                rl[R] = key[i];
                R++;
            }
            re[R - 1] = i;
            lastrun[key[i]] = R - 1;
        }
    if (!R)
        return 0;
    int cost = 0;
    /* stack of open pockets: class, gap sum, gap max */
    int sl[MAXL], ss[MAXL], sm[MAXL], sp = 0;
    int toplab[MAXL], topgap[MAXL], nt = 0;
    for (int j = 0; j < R; j++) {
        int X = rl[j];
        if (j > 0) {
            int g = rs[j] - re[j - 1] - 1 + (rs[j] == kink);
            if (sp) {
                ss[sp - 1] += g;
                if (g > sm[sp - 1])
                    sm[sp - 1] = g;
            } else
                topgap[nt - 1] = g;
        }
        if (sp && sl[sp - 1] == X) {
            cost += ss[sp - 1] - sm[sp - 1];
            sp--;
        } else if (!sp)
            toplab[nt++] = X;
        if (j < lastrun[X]) {
            sl[sp] = X;
            ss[sp] = sm[sp] = 0;
            sp++;
        }
    }
    /* top level */
    int cr[MAXL], m = 0;
    for (int t = 0; t < nt; t++)
        if (key[L + toplab[t] - 1] != PGHOST)
            cr[m++] = t;
    if (m) {
        for (int t = 0; t < cr[0]; t++)
            cost += topgap[t];
        for (int t = cr[m - 1]; t < nt - 1; t++)
            cost += topgap[t];
        int k = m / 2, best = 1 << 20;
        for (int i = 0; i < k; i++) {
            int sum = 0, mx = 0;
            /* left side: cr[i] .. cr[i+1] (or the middle) */
            int a = cr[i], b;
            if (i < k - 1)
                b = cr[i + 1];
            else
                b = cr[m - 1 - i];          /* middle: A_k .. (S ..) rho A_k */
            for (int t = a; t < b; t++) {
                sum += topgap[t];
                if (topgap[t] > mx)
                    mx = topgap[t];
            }
            if (i < k - 1)                  /* right side of ring i */
                for (int t = cr[m - 2 - i]; t < cr[m - 1 - i]; t++) {
                    sum += topgap[t];
                    if (topgap[t] > mx)
                        mx = topgap[t];
                }
            cost += sum - mx;
            if (mx < best)
                best = mx;
        }
        if (!(m & 1) && k)
            cost += best;                   /* no self class */
    } else {
        for (int t = 0; t < nt - 1; t++)
            cost += topgap[t];
    }
    if (!tch) {
        int d1 = rs[0], d2 = L - 1 - re[R - 1];
        cost += d1 < d2 ? d1 : d2;
    }
    return cost;
}

/* ------------------------------------------------------------------ */
static size_t MAXSTATES, MAXBYTES, MAXRECS, MAXKW;
static cnt_t *GF;            /* dense scratch counter vector */
static cnt_t *LT, *EQ;            /* results of one case, indexed by n */

static int uspan(int k) { return 2 * k + 1 - E; }

static void record(int ktop, const cnt_t *gf, int tch)
{
    int h = uspan(ktop);
    if (!tch || h < L)
        return;
    cnt_t *tgt = h == L ? EQ : LT;
    for (int wt = 0; wt <= NMAX; wt++)
        if (gf[wt])
            tgt[wt] = addm(tgt[wt], gf[wt]);
}

/* add gf (shifted by add) into table under key, pruned by the lower bound */
static void insert_pruned(table_t *t, const uint8_t *key, const cnt_t *gf,
                          int add, int lb)
{
    int lim = NMAX - lb - add;
    int any = 0;
    for (int wt = 0; wt <= lim; wt++)
        if (gf[wt]) {
            any = 1;
            break;
        }
    if (!any)
        return;
    size_t si = tab_get(t, key);
    for (int wt = 0; wt <= lim; wt++)
        if (gf[wt])
            tab_add(t, si, wt + add, gf[wt]);
}

/* one transition: decide position i of the current row */
static int step(work_t *w, const uint8_t *sk, int tch, int i, int occ, int last,
                int first_twin, uint8_t *key, int *ntp, int *mp)
{
    decode(w, sk);
    if (occ) {
        int nw = w->next++;
        w->partner[nw] = -2;
        w->f[i] = nw;
        for (int nb = i - 1; nb <= i + 1; nb += 2)
            if (nb >= 0 && nb < L && w->f[nb])
                unite(w, nw, w->f[nb]);
        if (first_twin) {
            for (int jj = L - 2 - i; jj <= L - i; jj += 2) {
                if (jj == i)
                    unite(w, nw, rho(w, nw));
                else if (jj >= 0 && jj < L && jj < i && w->f[jj])
                    unite(w, nw, rho(w, w->f[jj]));
            }
        }
    }
    int status = OK;
    if (i - 1 >= 0)
        status = retire(w, i - 1);
    if (status == OK && last && i + 1 < L)
        status = retire(w, i + 1);
    if (status != OK)
        return status;
    *ntp = tch || (occ && (i == 0 || i == L - 1));
    *mp = encode(w, key, *ntp);
    return OK;
}

/* lower bound (in weight) on what a nonempty state still needs; -1: never completes */
static int bound(const uint8_t *key, int nt, int k, int kmin, int i, int last, int first_twin)
{
    int lb = kmin > k ? 2 * (kmin - k) : 0;
    if (!nt) {
        int dist = L;
        for (int j = 0; j < L; j++)
            if (key[j]) {
                int dd = j < L - 1 - j ? j : L - 1 - j;
                if (dd < dist)
                    dist = dd;
            }
        if (2 * dist > lb)
            lb = 2 * dist;
    }
    if (!first_twin || last) {
        if (!completion_possible(key))
            return -1;
        int kink = (!last && !key[i] && i + 1 < L && key[i + 1]) ? i + 1 : -1;
        int eb = 2 * (OLDBOUND ? exact_bound(key, nt) : cactus_bound(key, nt, kink));
        if (eb > lb)
            lb = eb;
    }
    return lb;
}

static void count_case(int c, int d, int vhi)
{
    S = c + d;
    T = c - d;
    vlo = T - vhi;
    L = vhi - vlo + 1;
    E = S & 1;
    memset(LT, 0, NW * sizeof(cnt_t));
    memset(EQ, 0, NW * sizeof(cnt_t));
    if (L == 1) {
        if (c == 0 && d == 0)
            EQ[1] = addm(EQ[1], 1);
        return;
    }
    if (L > MAXL) {
        fprintf(stderr, "strip too wide\n");
        exit(1);
    }
    MAXC = (L + 1) / 2 + 1;
    GOFF = MAXC + 1;
    KL = L + MAXC + 1;
    KW = (3 * L + MAXC + 1 + 63) / 64;
    NW = NMAX + 1;

    int kmin = 0;
    while (uspan(kmin) < L)
        kmin++;

    table_t A, B, *cur = &A, *nxt = &B;
    tab_init(&A, 1024);
    tab_init(&B, 1024);
    uint8_t key[MAXL * 2 + 4];
    cnt_t one[MAXL * 2 + 4];
    work_t w;
    int u0;

    if (E == 0) {
        int um = S / 2;
        int half[MAXL], nh = 0;
        for (int i = 0; i < L; i++)
            if (((vlo + i - um) % 2 + 2) % 2 == 0 && 2 * (vlo + i) >= T)
                half[nh++] = i;
        /* all nonempty subsets of half with at most NMAX cells in the row */
        for (uint64_t mask = 1; mask < (1ull << nh); mask++) {
            int cells[MAXL], nc = 0, isc[MAXL];
            memset(isc, 0, sizeof(int) * L);
            for (int b = 0; b < nh; b++)
                if (mask >> b & 1) {
                    isc[half[b]] = 1;
                    isc[L - 1 - half[b]] = 1;
                }
            for (int i = 0; i < L; i++)
                if (isc[i])
                    cells[nc++] = i;
            if (nc > NMAX)
                continue;
            int lab[MAXL];
            for (int a = 0; a < MAXN; a++)
                w.parent[a] = a;
            for (int i = 0; i < L; i++)
                w.f[i] = 0;
            for (int a = 0; a < nc; a++) {
                lab[cells[a]] = a + 1;
                w.f[cells[a]] = a + 1;
            }
            for (int a = 0; a < nc; a++) {
                int i = cells[a], j = L - 1 - i;
                w.partner[a + 1] = j == i ? -1 : lab[j];
            }
            int tch = isc[0] || isc[L - 1];
            encode(&w, key, tch);
            memset(one, 0, sizeof(cnt_t) * NW);
            one[nc] = 1;
            insert_pruned(cur, key, one, 0, 0);
        }
        u0 = um + 1;
    } else {
        memset(key, 0, KL);
        memset(one, 0, sizeof(cnt_t) * NW);
        one[0] = 1;
        insert_pruned(cur, key, one, 0, 0);
        u0 = (S + 1) / 2;
    }

    int cells_u[MAXL];
    for (int k = 1; cur->n; k++) {
        int u = u0 + k - 1, nu = 0;
        for (int i = 0; i < L; i++)
            if (((vlo + i - u) % 2 + 2) % 2 == 0)
                cells_u[nu++] = i;
        int first_twin = (E == 1 && k == 1);
        for (int ci = 0; ci < nu; ci++) {
            int i = cells_u[ci];
            int last = ci == nu - 1;
            tab_clear(nxt);
            for (size_t si = 0; si < cur->n; si++) {
                uint8_t sk[MAXL * 2 + 4];
                unpack(cur->keys + si * KW, sk);
                const cnt_t *gf = GF;
                tab_load(cur, si, GF);
                int tch = sk[KL - 1];
                int empty = 1;
                for (int j = 0; j < L; j++)
                    if (sk[j]) {
                        empty = 0;
                        break;
                    }
                for (int occ = 0; occ <= 1; occ++) {
                    if (occ && k > 1 && empty)
                        continue;
                    int nt, m;
                    int status = step(&w, sk, tch, i, occ, last, first_twin, key, &nt, &m);
                    if (status == DEAD)
                        continue;
                    if (status == COMPLETE) {
                        record(k - 1, gf, tch);
                        continue;
                    }
                    int add = occ ? 2 : 0;
                    if (m == 0) {
                        if (first_twin && !last)
                            insert_pruned(nxt, key, gf, add, 0);
                        continue;
                    }
                    int lb = bound(key, nt, k, kmin, i, last, first_twin);
                    if (lb < 0)
                        continue;
                    if (LOOKAHEAD && !last && !first_twin && (LOOKAHEAD == 2 || (!key[i] && i + 1 < L && key[i + 1]))) {
                        /* one step ahead: the next cell is present or not */
                        int i2 = cells_u[ci + 1], last2 = ci + 1 == nu - 1;
                        int la = 1 << 20;
                        uint8_t k2[MAXL * 2 + 4];
                        for (int o2 = 0; o2 <= 1 && la; o2++) {
                            int nt2, m2;
                            int st2 = step(&w, key, nt, i2, o2, last2, 0, k2, &nt2, &m2);
                            if (st2 == DEAD)
                                continue;
                            int val = 2 * o2;
                            if (st2 == OK) {
                                if (m2 == 0)
                                    continue;
                                int lb2 = bound(k2, nt2, k, kmin, i2, last2, 0);
                                if (lb2 < 0)
                                    continue;
                                val += lb2;
                            }
                            if (val < la)
                                la = val;
                        }
                        if (la > lb)
                            lb = la;
                    }
                    insert_pruned(nxt, key, gf, add, lb);
                }
                tab_unload(cur, si, GF);
            }
            table_t *x = cur;
            cur = nxt;
            nxt = x;
            if (cur->n > MAXSTATES)
                MAXSTATES = cur->n;
            size_t by = tab_bytes(cur) + tab_bytes(nxt);
            if (by > MAXBYTES)
                MAXBYTES = by, MAXRECS = cur->nrec, MAXKW = KW;
        }
        if (REFLECT && T == 0 && cur->n) {
            /* merge each state with its mirror image v -> -v */
            tab_clear(nxt);
            uint8_t rkey[MAXL * 2 + 4];
            for (size_t si = 0; si < cur->n; si++) {
                uint8_t sk[MAXL * 2 + 4];
                unpack(cur->keys + si * KW, sk);
                decode(&w, sk);
                for (int a = 0, b = L - 1; a < b; a++, b--) {
                    int x = w.f[a];
                    w.f[a] = w.f[b];
                    w.f[b] = x;
                }
                encode(&w, rkey, sk[KL - 1]);
                const uint8_t *use = memcmp(rkey, sk, KL) < 0 ? rkey : sk;
                tab_load(cur, si, GF);
                insert_pruned(nxt, use, GF, 0, 0);
                tab_unload(cur, si, GF);
            }
            table_t *x = cur;
            cur = nxt;
            nxt = x;
        }
    }
    tab_free(&A);
    tab_free(&B);
}

int main(int argc, char **argv)
{
    int a = 1;
    while (a < argc && argv[a][0] == '-') {
        if (!strcmp(argv[a], "-m") && a + 1 < argc)
            MOD = strtoull(argv[++a], 0, 10), a++;
        else if (!strcmp(argv[a], "-A"))
            LOOKAHEAD = 2, a++;
        else if (!strcmp(argv[a], "-N"))
            LOOKAHEAD = 0, a++;
        else if (!strcmp(argv[a], "-K"))
            OLDBOUND = 1, a++;
        else if (!strcmp(argv[a], "-C"))
            CHECK = 1, a++;
        else if (!strcmp(argv[a], "-R"))
            REFLECT = 0, a++;
        else
            break;
    }
    if (a >= argc) {
        fprintf(stderr, "usage: %s [-m M] [-R] N [c d vhi]\n", argv[0]);
        return 1;
    }
    NMAX = atoi(argv[a]);
    NW = NMAX + 1;
    LT = calloc(NW, sizeof(cnt_t));
    GF = calloc(NW, sizeof(cnt_t));
    EQ = calloc(NW, sizeof(cnt_t));
    clock_t t0 = clock();
    if (argc - a == 4) {
        int c = atoi(argv[a + 1]), d = atoi(argv[a + 2]), vhi = atoi(argv[a + 3]);
        count_case(c, d, vhi);
        for (int n = 1; n <= NMAX; n++)
            printf("%d %llu %llu\n", n, (unsigned long long)LT[n],
                   (unsigned long long)EQ[n]);
    } else {
        cnt_t *out = calloc(NW, sizeof(cnt_t));
        for (int c = 0; c <= 1; c++)
            for (int d = 0; d <= 1; d++) {
                int t = c - d;
                int vhi = t >= 0 ? (t + 1) / 2 : 0;
                for (; 2 * vhi - t + 1 <= NMAX; vhi++) {
                    if (2 * vhi - t + 1 < 1)
                        continue;
                    count_case(c, d, vhi);
                    for (int n = 1; n <= NMAX; n++) {
                        out[n] = addm(out[n], addm(LT[n], LT[n]));
                        out[n] = addm(out[n], EQ[n]);
                    }
                }
            }
        for (int n = 1; n <= NMAX; n++)
            printf("%d %llu\n", n, (unsigned long long)out[n]);
    }
    fprintf(stderr, "# peak memory %.1f MB (both tables; %zu count records, key %zu words)\n",
            MAXBYTES / 1e6, MAXRECS, MAXKW);
    fprintf(stderr, "# max states %zu, %.2fs\n", MAXSTATES,
            (double)(clock() - t0) / CLOCKS_PER_SEC);
    return 0;
}
