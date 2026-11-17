This is a CWEB change file for diag180.w.  It produces the production
version of the program, without the experimental options -K (old bound),
-N and -B (variants of the lookahead), -H (no height terms), and -C
(check of packing and unpacking).

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
@ Two states are compared lexicographically.  (If the labels agree, the
number of components is the same.)  This is needed only for the check below.
@<Sub...@>=
static int compare_states(const state_t *a, const state_t *b)
{
    /* on entry: both states are canonical; only |label[0..W-1]| and |partner[1..num_components]| are defined */
    int c = memcmp(a->label, b->label, W);
    if (c)
        return c;
    c = memcmp(a->partner + 1, b->partner + 1, a->num_components);
    if (c)
        return c;
    return a->touches - b->touches;
}

@ The function |pack_and_check| packs a state.  With the option \texttt{-C},
it checks that unpacking gives the state back (we have run this check for $N\le 30$).
@<Sub...@>=
static void pack_and_check(const state_t *s, int decided_position, uint8_t *packed_key)
{
    pack(s, decided_position, packed_key);
    if (CHECK_PACKING) {
        state_t unpacked;
        unpack(packed_key, decided_position, &unpacked);
        if (unpacked.num_components != s->num_components || compare_states(&unpacked, s)) {
            fprintf(stderr, "pack/unpack mismatch\n");
            exit(2);
        }
    }
}
@y
@ The function |pack_and_check| packs a state.  (In the full version of the program,
the option \texttt{-C} checks that unpacking gives the state back.)
@<Sub...@>=
static void pack_and_check(const state_t *s, int decided_position, uint8_t *packed_key)
{
    pack(s, decided_position, packed_key);
}
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
@<Set |connection_term| to the old spanning tree bound@>=
        int representative[MAXW + 1];
        int occupied[MAXW], num_occupied = 0;
        for_int_from_to(j, 1, state->num_components) {
            int q = state->partner[j];
            representative[j] = (q != 0 && q < j) ? q : j;
        }
        for_int_from_to(j, 0, W - 1)
            if (label[j])
                occupied[num_occupied++] = j;
        /* Kruskal over gaps between consecutive occupied positions */
        int parent[MAXW + 1];
        int num_sets = 0;
        int seen[MAXW + 1];
        memset(seen, 0, sizeof(seen));
        for_int_from_to(a, 0, num_occupied - 1) {
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
            for_int_from_to(a, 0, num_gaps - 1) {
                gap[a] = occupied[a + 1] - occupied[a] - 1;
                bucket[gap[a]]++;
            }
            for_int_from_to(g, 1, W)
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
        connection_term = 2 * cost;

@y
@z

@x
    int connection_term;
    if (USE_OLD_BOUND) {
        @<Set |connection_term| to the old spanning tree bound@>@;
    } else if (state->touches && frontier_height >= TARGET_ROW)
        connection_term = 2 * connection_bound(state);
    else
        connection_term = 2 * cactus_bound(state);
@y
    int connection_term;
    if (state->touches && frontier_height >= TARGET_ROW)
        connection_term = 2 * connection_bound(state);
    else
        connection_term = 2 * cactus_bound(state);
@z

@x
(With the option \texttt{-N} the lookahead is switched off, and with \texttt{-B}
the old version of the lookahead is used, which tries both possibilities for the
kink cell with |step| and takes the smaller of the two bounds.)
@y
@z

@x
                        if (LOOKAHEAD && BLOCKED_LOOKAHEAD && !last_in_row &&
@y
                        if (!last_in_row &&
@z

@x
                        } else {
                            lowerbound = bound(&new_state, k, min_rows, i, last_in_row);
                            if (lowerbound >= 0 && LOOKAHEAD && !last_in_row &&
                                !new_state.label[i] && new_state.label[i + 1]) {
                                @<Sharpen |lowerbound| by trying both possibilities for the next cell@>@;
                            }
                        }
@y
                        } else
                            lowerbound = bound(&new_state, k, min_rows, i, last_in_row);
@z

@x
@ {\bf Trying both possibilities.}
With the option \texttt{-B}, both possibilities for the kink cell are tried with
|step|, and the bound is computed for both resulting states.
A completion counts only if it would be a valid polyomino.
@<Sharpen |lowerbound| by trying both...@>=
int next_i = row_cells[cell_number + 1], next_is_last = cell_number + 1 == num_row_cells - 1;
int lookahead_bound = 1 << 20;
state_t state_ahead;
for (int next_occupied = 0; next_occupied <= 1 && lookahead_bound; next_occupied++) {
    int status_ahead = step(&new_state, next_i, next_occupied, next_is_last, &state_ahead);
    if (status_ahead == DEAD)
        continue;
    /* a completion only counts if it is a valid polyomino */
    if (status_ahead == COMPLETE && !(new_state.touches && vspan(k - 1) >= W))
        continue;
    int bound_this_way = 2 * next_occupied;
    if (status_ahead == OK) {
        int lowerbound_ahead = bound(&state_ahead, k, min_rows, next_i, next_is_last);
        if (lowerbound_ahead < 0)
            continue;
        bound_this_way += lowerbound_ahead;
    }
    lookahead_bound = min(lookahead_bound, bound_this_way);
}
lowerbound = max(lowerbound, lookahead_bound);

@y
@z

@x
  -N     no lookahead
  -H     no height terms in the bound
  -K     use the old bound (minimum spanning tree only)
  -C     check that packing and unpacking of every key are consistent
  -B     at the kink, try both choices instead of blocking the gaps
@y
@z

@x
static int USE_OLD_BOUND = 0; /* \texttt{-K} */
@y
@z

@x
static int LOOKAHEAD = 1; /* lookahead at the kink? (option \texttt{-N} turns it off) */
static int BLOCKED_LOOKAHEAD = 1; /* use blocked gaps instead of trying both choices (option \texttt{-B} turns it off) */
static int CHECK_PACKING = 0; /* \texttt{-C} */
@y
@z

@x
        else if (!strcmp(argv[arg], "-N"))
            LOOKAHEAD = 0, arg++;
        else if (!strcmp(argv[arg], "-H"))
            USE_HEIGHT_TERMS = 0, arg++;
        else if (!strcmp(argv[arg], "-K"))
            USE_OLD_BOUND = 1, arg++;
        else if (!strcmp(argv[arg], "-C"))
            CHECK_PACKING = 1, arg++;
        else if (!strcmp(argv[arg], "-B"))
            BLOCKED_LOOKAHEAD = 0, arg++;
@y
@z

@x
        fprintf(stream, "# bound: %s, height term %s, lookahead %s; reflection merge %s; pack check %s\n",
                USE_OLD_BOUND ? "old Kruskal (-K)" : "cactus", USE_HEIGHT_TERMS ? "on" : "off (-H)",
                !LOOKAHEAD ? "off (-N)" : BLOCKED_LOOKAHEAD ? "at the kink, blocked gaps" : "at the kink, both choices (-B)",
                MERGE_MIRROR_STATES ? "on" : "off (-R)", CHECK_PACKING ? "on (-C)" : "off");
@y
        fprintf(stream, "# production version (cactus bound with height terms, lookahead at the kink); reflection merge %s\n",
                MERGE_MIRROR_STATES ? "on" : "off (-R)");
@z

@x
separately.  Different variants of the bound (\texttt{-K}, \texttt{-H},
\texttt{-N}, \texttt{-B}) lead to the same counts; this is a strong test, since an
invalid bound would lose polyominoes.
@y
separately.  In the full version of the program, different variants of
the bound and of the lookahead can be selected, and they lead to the same
counts; this is a strong test, since an
invalid bound would lose polyominoes.
@z
