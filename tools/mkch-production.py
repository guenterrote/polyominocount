# one-time generator for src/production.ch (exact copies of master lines)
L=open('src/diag180.w').read().split('\n')
def find(prefix, start=0):
    for i in range(start,len(L)):
        if L[i].startswith(prefix): return i
    raise SystemExit("not found: "+prefix)
out=["""This is a CWEB change file for diag180.w.  It produces the production
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
i=find("static int USE_HEIGHT_TERMS = 1;"); ch(i,i+1,"")
i=find("        return USE_HEIGHT_TERMS ? missing_rows : 0;"); assert L[i+1]=="    if (!USE_HEIGHT_TERMS)" and L[i+2]=="        return SIDE_DISTANCE;"
ch(i,i+3,"        return missing_rows;")
i=find("@ An older and weaker bound is retained"); j=find("@ The function |bound| returns"); ch(i,j,"")
i=find("    int connection_bound = 2 * (USE_OLD_BOUND ?"); ch(i,i+1,"    int connection_bound = 2 * cactus_bound(state, kink);")
i=find("the lookahead counts only if it would be a valid polyomino.  With the")
assert L[i+2]=="is switched off." and L[i+3].startswith("@<Sharpen") and L[i+4].lstrip().startswith("if (LOOKAHEAD &&")
ch(i,i+5,"the lookahead counts only if it would be a valid polyomino.\n"+L[i+3]+"\n"
   "                    if (!last_in_row && !new_state.label[i] && i + 1 < W && new_state.label[i + 1]) {")
i=find("  -N     no lookahead;"); assert L[i+2].startswith("  -K")
ch(i,i+3,"")
i=find("static int USE_OLD_BOUND = 0;"); ch(i,i+1,"")
i=find("static int LOOKAHEAD = 1;"); ch(i,i+1,"")
i=find('        else if (!strcmp(argv[arg], "-A"))'); assert L[i+7]=="            USE_OLD_BOUND = 1, arg++;"
ch(i,i+8,"")
i=find('        fprintf(stream, "# bound: %s, height term'); assert "MERGE_MIRROR_STATES" in L[i+3]
ch(i,i+4,'        fprintf(stream, "# production version (cactus bound with height terms, lookahead at the kink); reflection merge %s; pack check %s\\n",\n'
   '                MERGE_MIRROR_STATES ? "on" : "off (-R)", CHECK_PACKING ? "on (-C)" : "off");')
i=find("separately.  Different variants of the bound"); assert L[i+2].startswith("invalid bound would lose")
ch(i,i+3,"separately.  In the full version of the program, different variants of\nthe bound and of the lookahead can be selected, and they lead to the same\ncounts; this is a strong test, since an\ninvalid bound would lose polyominoes.")
open('src/production.ch','w').write("\n".join(out))
