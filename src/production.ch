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
static int USE_HEIGHT_TERMS = 1; /* use the height terms? (option \texttt{-H} turns them off) */
@y
@z

@x
        return USE_HEIGHT_TERMS ? missing_rows : 0;
    if (!USE_HEIGHT_TERMS)
        return SIDE_DISTANCE;
@y
        return missing_rows;
@z

@x
@*1 The old bound.
An older and weaker bound is retained for comparison (option
\texttt{-K}).  It computes a minimum spanning tree over the gaps by
Kruskal's algorithm, identifying partners, and adds the distance to the
side.  It ignores the kink, the central term, and the height.
The code is part of the function |bound| below; the state is not empty.
@<Set |connection_bound| to the old spanning tree bound@>=
        int representative[MAXW + 1];
        int occupied[MAXW], num_occupied = 0;
        for (int j = 1; j <= state->num_components; j++) {
            int q = state->partner[j];
            representative[j] = (q != 0 && q < j) ? q : j;
        }
        for (int j = 0; j < W; j++)
            if (label[j])
                occupied[num_occupied++] = j;
        /* Kruskal over gaps between consecutive occupied positions */
        int parent[MAXW + 1];
        int num_sets = 0;
        int seen[MAXW + 1];
        memset(seen, 0, sizeof(seen));
        for (int a = 0; a < num_occupied; a++) {
            int rep = representative[label[occupied[a]]];
            if (!seen[rep]) {
                seen[rep] = 1;
                num_sets++;
            }
            parent[rep] = rep;
        }
        int cost = 0;
        if (num_sets > 1) {
            /* counting sort of gaps by weight */
            int num_gaps = num_occupied - 1;
            int gap[MAXW], order[MAXW], bucket[MAXW + 1];
            memset(bucket, 0, sizeof(bucket));
            for (int a = 0; a < num_gaps; a++) {
                gap[a] = occupied[a + 1] - occupied[a] - 1;
                bucket[gap[a]]++;
            }
            for (int g = 1; g <= W; g++)
                bucket[g] += bucket[g - 1];
            for (int a = num_gaps - 1; a >= 0; a--)
                order[--bucket[gap[a]]] = a;
            for (int e = 0; e < num_gaps && num_sets > 1; e++) {
                int a = order[e];
                int x = representative[label[occupied[a]]], y = representative[label[occupied[a + 1]]];
                while (parent[x] != x)
                    x = parent[x];
                while (parent[y] != y)
                    y = parent[y];
                if (x != y) {
                    parent[x] = y;
                    cost += gap[a];
                    num_sets--;
                }
            }
        }
        if (!state->touches) {
            int left_distance = occupied[0], right_distance = W - 1 - occupied[num_occupied - 1];
            cost += min(left_distance, right_distance);
        }
        connection_bound = 2 * cost;

@y
@z

@x
    int connection_bound;
    if (USE_OLD_BOUND) {
        @<Set |connection_bound| to the old spanning tree bound@>@;
    } else
        connection_bound = 2 * cactus_bound(state, kink);
@y
    int connection_bound = 2 * cactus_bound(state, kink);
@z

@x
the lookahead counts only if it would be a valid polyomino.  With the
option \texttt{-A} the lookahead is done for every cell, and with \texttt{-N} it
is switched off.
@<Sharpen |lowerbound|...@>=
                    if (LOOKAHEAD && !last_in_row && (LOOKAHEAD == 2 || (!new_state.label[i] && i + 1 < W && new_state.label[i + 1]))) {
@y
the lookahead counts only if it would be a valid polyomino.
@<Sharpen |lowerbound|...@>=
                    if (!last_in_row && !new_state.label[i] && i + 1 < W && new_state.label[i + 1]) {
@z

@x
  -N     no lookahead;   -A  lookahead after every cell
  -H     no height terms in the bound
  -K     use the old bound (minimum spanning tree only)
@y
@z

@x
static int USE_OLD_BOUND = 0; /* \texttt{-K} */
@y
@z

@x
static int LOOKAHEAD = 1; /* 0 none, 1 at the kink, 2 after every cell */
@y
@z

@x
        else if (!strcmp(argv[arg], "-A"))
            LOOKAHEAD = 2, arg++;
        else if (!strcmp(argv[arg], "-N"))
            LOOKAHEAD = 0, arg++;
        else if (!strcmp(argv[arg], "-H"))
            USE_HEIGHT_TERMS = 0, arg++;
        else if (!strcmp(argv[arg], "-K"))
            USE_OLD_BOUND = 1, arg++;
@y
@z

@x
        fprintf(stream, "# bound: %s, height term %s, lookahead %s; reflection merge %s; pack check %s\n",
                USE_OLD_BOUND ? "old Kruskal (-K)" : "cactus", USE_HEIGHT_TERMS ? "on" : "off (-H)",
                LOOKAHEAD == 2 ? "every cell (-A)" : LOOKAHEAD ? "at the kink" : "off (-N)",
                MERGE_MIRROR_STATES ? "on" : "off (-R)", CHECK_PACKING ? "on (-C)" : "off");
@y
        fprintf(stream, "# production version (cactus bound with height terms, lookahead at the kink); reflection merge %s; pack check %s\n",
                MERGE_MIRROR_STATES ? "on" : "off (-R)", CHECK_PACKING ? "on (-C)" : "off");
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
