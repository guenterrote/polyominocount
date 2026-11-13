This is a CWEB change file for diag180.w.  It produces the production
version of the program, without the experimental options -K (old bound),
-A and -N (variants of the lookahead), and -H (no height terms).

    ctangle diag180.w production.ch diag180p.c      (make tangle)
    cweave  diag180.w production.ch diag180p.tex    (make doc)

Everything else, in particular the explanations, comes from diag180.w.
Each change below replaces the lines between @x and @y, which must match
the lines of diag180.w exactly, by the lines between @y and @z.  If one
of these places is edited in diag180.w, ctangle and cweave report
"Change file entry did not match", and the lines after @x must be
updated here.

@x
\title{\textbf{diag180}: Counting polyominoes with $180^\circ$ rotational
  symmetry\\ by a transfer matrix along diagonals}
@y
\title{\textbf{diag180}: Counting polyominoes with $180^\circ$ rotational
  symmetry\\ by a transfer matrix along diagonals\\ \large (production version)}
@z

@x
static int HEIGHT = 1; /* use the height terms? (option \texttt{-H} turns them off) */
@y
@z

@x
        return HEIGHT ? top : 0;
    if (!HEIGHT)
        return DEDGE;
@y
        return top;
@z

@x
@ An older and weaker bound is retained for comparison (option
\texttt{-K}).  It computes a minimum spanning tree over the gaps by
Kruskal's algorithm, identifying partners, and adds the distance to the
side.  It ignores the kink, the central term, and the height.
@<Sub...@>=
static int exact_bound(const state_t *s)
{
    const uint8_t *label = s->label;
    int touches = s->touches;
    int sup[MAXW + 1];
    int occ[MAXW], no = 0;
    for (int j = 1; j <= s->num_components; j++) {
        int q = s->partner[j];
        sup[j] = (q != 0 && q < j) ? q : j;
    }
    for (int i = 0; i < W; i++)
        if (label[i])
            occ[no++] = i;
    if (!no)
        return 0;
    /* Kruskal over gaps between consecutive occupied positions */
    int par[MAXW + 1];
    int ncls = 0;
    int seen[MAXW + 1];
    memset(seen, 0, sizeof(seen));
    for (int a = 0; a < no; a++) {
        int s = sup[label[occ[a]]];
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
        int wgt[MAXW], ord[MAXW], bucket[MAXW + 1];
        memset(bucket, 0, sizeof(bucket));
        for (int a = 0; a < ne; a++) {
            wgt[a] = occ[a + 1] - occ[a] - 1;
            bucket[wgt[a]]++;
        }
        for (int g = 1; g <= W; g++)
            bucket[g] += bucket[g - 1];
        for (int a = ne - 1; a >= 0; a--)
            ord[--bucket[wgt[a]]] = a;
        for (int e = 0; e < ne && ncls > 1; e++) {
            int a = ord[e];
            int x = sup[label[occ[a]]], y = sup[label[occ[a + 1]]];
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
    if (!touches) {
        int d1 = occ[0], d2 = W - 1 - occ[no - 1];
        cost += min(d1, d2);
    }
    return cost;
}

@y
@z

@x
    int eb = 2 * (OLDBOUND ? exact_bound(s) : cactus_bound(s, kink));
@y
    int eb = 2 * cactus_bound(s, kink);
@z

@x
the lookahead counts only if it would be a valid polyomino.  With the
option \texttt{-A} the lookahead is done for every cell, and with \texttt{-N} it
is switched off.
@<Sharpen |lowerbound|...@>=
                    if (LOOKAHEAD && !last && (LOOKAHEAD == 2 || (!key.label[i] && i + 1 < W && key.label[i + 1]))) {
@y
the lookahead counts only if it would be a valid polyomino.
@<Sharpen |lowerbound|...@>=
                    if (!last && !key.label[i] && i + 1 < W && key.label[i + 1]) {
@z

@x
  -N     no lookahead;   -A  lookahead after every cell
  -H     no height terms in the bound
  -K     use the old bound (minimum spanning tree only)
@y
@z

@x
static int OLDBOUND = 0; /* \texttt{-K} */
@y
@z

@x
static int LOOKAHEAD = 1; /* 0 none, 1 at the kink, 2 after every cell */
@y
@z

@x
        else if (!strcmp(argv[a], "-A"))
            LOOKAHEAD = 2, a++;
        else if (!strcmp(argv[a], "-N"))
            LOOKAHEAD = 0, a++;
        else if (!strcmp(argv[a], "-H"))
            HEIGHT = 0, a++;
        else if (!strcmp(argv[a], "-K"))
            OLDBOUND = 1, a++;
@y
@z

@x
        fprintf(f, "# bound: %s, height term %s, lookahead %s; reflection merge %s; pack check %s\n",
                OLDBOUND ? "old Kruskal (-K)" : "cactus", HEIGHT ? "on" : "off (-H)",
                LOOKAHEAD == 2 ? "every cell (-A)" : LOOKAHEAD ? "at the kink" : "off (-N)",
                REFLECT ? "on" : "off (-R)", CHECK ? "on (-C)" : "off");
@y
        fprintf(f, "# production version (cactus bound with height terms, lookahead at the kink); reflection merge %s; pack check %s\n",
                REFLECT ? "on" : "off (-R)", CHECK ? "on (-C)" : "off");
@z

@x
separately.  Different variants of the bound (\texttt{-K}, \texttt{-H},
\texttt{-N}, \texttt{-A}) lead to the same counts; this is a strong test, since an
invalid bound would lose polyominoes.
@y
separately.  In the full version of the program, different variants of
the bound and of the lookahead can be selected, and they lead to the same
counts; this is a strong test, since an
invalid bound would lose polyominoes.
@z
