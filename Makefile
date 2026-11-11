CC = gcc
CFLAGS = -O3 -march=native -Wall
# make CNT128=1: 128-bit counters in diag180 (exact, no CRT needed up to n = 71)
ifdef CNT128
D180FLAGS = -DCNT128
endif

all: build/redelmeier build/ustats build/diag180 build/diag180p

build/ustats: src/ustats.c
	mkdir -p build
	$(CC) $(CFLAGS) -o $@ $<

build/redelmeier: src/redelmeier.c
	mkdir -p build
	$(CC) $(CFLAGS) -o $@ $<

build/diag180: src/diag180.c FORCE180
	mkdir -p build
	$(CC) $(CFLAGS) $(D180FLAGS) -o $@ $<

build/diag180p: src/diag180p.c FORCE180
	mkdir -p build
	$(CC) $(CFLAGS) $(D180FLAGS) -o $@ $<

test: build/redelmeier build/diag180 build/diag180p
	sh tests/run_tests.sh
	sh tests/test_mirror.sh
	sh tests/test_diag180.sh

# diag180 is a CWEB program: src/diag180.w is the source.
# src/production.ch is a CWEB change file that turns it into the production
# version diag180p (without the experimental options -K, -A, -N, -H).
# "make tangle" regenerates src/diag180.c and src/diag180p.c (which are kept
# in the repository, so that ctangle is not needed for compiling);
# "make doc" produces doc/diag180.pdf and doc/diag180p.pdf (needs cweave,
# pdflatex and bibtex with the style plainurl; the bibliography is
# src/diag180.bib; the LaTeX class cweb-hy is in cweb/).
tangle:
	cd src && ctangle diag180.w - diag180.c
	cd src && ctangle diag180.w production.ch diag180p.c

LATEX = TEXINPUTS=../cweb: pdflatex -interaction=nonstopmode -halt-on-error

doc:
	mkdir -p doc
	cd doc && cweave ../src/diag180.w - diag180.tex
	cd doc && cweave ../src/diag180.w ../src/production.ch diag180p.tex
	for d in diag180 diag180p; do ( cd doc && \
	  $(LATEX) $$d.tex >/dev/null && \
	  BIBINPUTS=../src: bibtex $$d >/dev/null && \
	  $(LATEX) $$d.tex >/dev/null && \
	  $(LATEX) $$d.tex >/dev/null && \
	  $(LATEX) $$d.tex | grep -E "^Output|Warning.*undefined" ; \
	  rm -f $$d.aux $$d.log $$d.out $$d.toc $$d.idx $$d.scn $$d.tex $$d.bbl $$d.blg ) ; done

# always rebuild diag180, so switching CNT128 on or off takes effect
FORCE180:

.PHONY: test FORCE180 tangle doc
