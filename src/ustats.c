/*
 * ustats.c -- statistics of Redelmeier's untried set near the leaves.
 *
 * Usage: ustats N K
 * For every node of Redelmeier's recursion (stack discipline) where the
 * current polyomino S has N-K cells (K cells still missing), record
 *   |U|      size of the untried set available to the subtree
 *   L1..L3   number of "unseen" cells (not in S, U, X or the forbidden
 *            half plane) at BFS distance 1, 2, 3 from U, moving only
 *            through unseen cells.
 * Prints the distribution of |U| and the averages of L1..L3.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 40
typedef unsigned long long u64;

static int N, K, W;
static unsigned char *mark;
static int *untried[MAXN + 1];
static int off[4];

static u64 hist[4 * MAXN + 8], nodes;
static double sumL[4];
static int *dist, *bfsq;

static void stats(int *u, int nu)
{
    nodes++;
    hist[nu]++;
    int qh = 0, qt = 0;
    for (int i = 0; i < nu; i++) { dist[u[i]] = 0; bfsq[qt++] = u[i]; }
    int cnt[4] = {0, 0, 0, 0};
    while (qh < qt) {
        int c = bfsq[qh++];
        if (dist[c] >= 3) continue;
        for (int j = 0; j < 4; j++) {
            int d = c + off[j];
            if (!mark[d] && dist[d] < 0) { dist[d] = dist[c] + 1; cnt[dist[d]]++; bfsq[qt++] = d; }
        }
    }
    for (int i = 0; i < qt; i++) dist[bfsq[i]] = -1;
    for (int l = 1; l <= 3; l++) sumL[l] += cnt[l];
}

static void rec(int depth, int *u, int nu)
{
    if (depth == N - K) { stats(u, nu); return; }
    int *ns = untried[depth + 1];
    while (nu > 0) {
        int c = u[--nu];
        int k = depth + 1;
        memcpy(ns, u, nu * sizeof(int));
        int m = nu;
        for (int j = 0; j < 4; j++) {
            int d = c + off[j];
            if (!mark[d]) { mark[d] = 1; ns[m++] = d; }
        }
        rec(k, ns, m);
        for (int j = nu; j < m; j++) mark[ns[j]] = 0;
    }
}

int main(int argc, char **argv)
{
    if (argc < 3) { fprintf(stderr, "usage: %s N K\n", argv[0]); return 1; }
    N = atoi(argv[1]); K = atoi(argv[2]);
    W = 2 * N + 8;
    int H = N + 8;
    mark = calloc((size_t)W * H, 1);
    dist = malloc(sizeof(int) * W * H);
    bfsq = malloc(sizeof(int) * W * H);
    for (int i = 0; i < W * H; i++) dist[i] = -1;
    for (int i = 0; i <= N; i++) untried[i] = malloc(sizeof(int) * (4 * N + 8));
    off[0] = 1; off[1] = -1; off[2] = W; off[3] = -W;
    /* guard rows/columns so BFS never leaves the array */
    for (int y = 0; y < H; y++) for (int x = 0; x < W; x++)
        if (y < 4 || y >= H - 2 || x < 2 || x >= W - 2) mark[y * W + x] = 1;
    int oy = 4, ox = W / 2;
    for (int x = 0; x < ox; x++) mark[oy * W + x] = 1;   /* y = 0, x < 0 forbidden */
    int origin = oy * W + ox;
    mark[origin] = 1;
    untried[0][0] = origin;
    rec(0, untried[0], 1);
    printf("N=%d, %d missing cells: %llu nodes with |S|=%d\n", N, K, nodes, N - K);
    double s = 0;
    for (int i = 0; i < 4 * MAXN + 8; i++) s += (double)i * hist[i];
    printf("average |U| = %.2f\n", s / nodes);
    printf("average layer sizes at distance 1,2,3 from U: %.2f %.2f %.2f\n",
           sumL[1] / nodes, sumL[2] / nodes, sumL[3] / nodes);
    printf("|U| histogram:");
    for (int i = 0; i < 4 * MAXN + 8; i++) if (hist[i]) printf(" %d:%.1f%%", i, 100.0 * hist[i] / nodes);
    printf("\n");
    return 0;
}
