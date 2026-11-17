# one-time generator for src/production.ch (exact copies of master lines)
L=open('src/diag180.w').read().split('\n')
def find(prefix, start=0):
    for i in range(start,len(L)):
        if L[i].startswith(prefix): return i
    raise SystemExit("not found: "+prefix)
out=["""This is a CWEB change file for diag180.w.  It produces the production
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
"""]
pos=0
def ch(i,j,new):
    global pos
    assert i>=pos, (i,pos)
    pos=j
    out.append("@x\n"+"\n".join(L[i:j])+"\n@y\n"+(new+"\n" if new else "")+"@z\n")
i=find(r"\title{\textbf{diag180}"); 
ch(i,i+2, L[i]+"\n"+L[i+1].replace("diagonals}","diagonals\\\\ \\large (production version)}"))
assert "diagonals}" in L[i+1]
h=find("@ Two states are compared lexicographically.")
i=find("@ The function |pack_and_check| packs a state.  With the option"); assert L[i+5]=="    pack(s, decided_position, packed_key);" and L[i+14]=="}"
ch(h,i+15,"@ The function |pack_and_check| packs a state.  (In the full version of the program,\n"
   "the option \\texttt{-C} checks that unpacking gives the state back.)\n"+"\n".join(L[i+2:i+6])+"\n}")
i=find("static int USE_HEIGHT_TERMS = 1;"); ch(i,i+1,"")
i=find("        return USE_HEIGHT_TERMS ? missing_rows : 0;"); assert L[i+1]=="    if (!USE_HEIGHT_TERMS)" and L[i+2]=="        return SIDE_DISTANCE;"
ch(i,i+3,"        return missing_rows;")
i=find("@*1 The old bound."); j=find("@*1 The complete bound."); ch(i,j,"")
i=find("    int connection_term;"); assert L[i+6]=="        connection_term = 2 * cactus_bound(state);"
ch(i,i+7,"    int connection_term;\n    if (state->touches && frontier_height >= TARGET_ROW)\n"
   "        connection_term = 2 * connection_bound(state);\n    else\n        connection_term = 2 * cactus_bound(state);")
i=find("the lookahead counts only if it would be a valid polyomino.  With the")
assert L[i+1]=="option \\texttt{-N} the lookahead is switched off."
ch(i,i+2,"the lookahead counts only if it would be a valid polyomino.")
i=find("                        if (LOOKAHEAD && cell_number + 1 < num_row_cells - 1 &&")
ch(i,i+1,"                        if (cell_number + 1 < num_row_cells - 1 &&")
i=find("                        } else if (kink_cell_optional && BLOCKED_LOOKAHEAD) {")
ch(i,i+1,"                        } else if (kink_cell_optional) {")
i=find("                            if (lowerbound >= 0 && LOOKAHEAD && !last_in_row &&")
assert L[i+1].lstrip().startswith("!new_state.label[i] && i + 1 < W")
ch(i,i+2,"                            if (lowerbound >= 0 && !last_in_row && !new_state.label[i] && new_state.label[i + 1]) {")
i=find("  -N     no lookahead"); assert L[i+2].startswith("  -K") and L[i+3].startswith("  -C") and L[i+4].startswith("  -B")
ch(i,i+5,"")
i=find("static int USE_OLD_BOUND = 0;"); ch(i,i+1,"")
i=find("static int LOOKAHEAD = 1;"); assert L[i+1].startswith("static int BLOCKED_LOOKAHEAD") and L[i+2].startswith("static int CHECK_PACKING"); ch(i,i+3,"")
i=find('        else if (!strcmp(argv[arg], "-N"))'); assert L[i+5]=="            USE_OLD_BOUND = 1, arg++;" and L[i+7]=="            CHECK_PACKING = 1, arg++;" and L[i+9]=="            BLOCKED_LOOKAHEAD = 0, arg++;"
ch(i,i+10,"")
i=find('        fprintf(stream, "# bound: %s, height term'); assert "MERGE_MIRROR_STATES" in L[i+3]
ch(i,i+4,'        fprintf(stream, "# production version (cactus bound with height terms, lookahead at the kink); reflection merge %s\\n",\n'
   '                MERGE_MIRROR_STATES ? "on" : "off (-R)");')
i=find("separately.  Different variants of the bound"); assert L[i+2].startswith("invalid bound would lose")
ch(i,i+3,"separately.  In the full version of the program, different variants of\nthe bound and of the lookahead can be selected, and they lead to the same\ncounts; this is a strong test, since an\ninvalid bound would lose polyominoes.")
open('src/production.ch','w').write("\n".join(out))
