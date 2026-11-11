CC = gcc
CFLAGS = -O3 -march=native -Wall
# make CNT128=1: 128-bit counters in diag180 (exact, no CRT needed up to n = 71)
ifdef CNT128
D180FLAGS = -DCNT128
endif

all: build/redelmeier build/ustats build/diag180

build/ustats: src/ustats.c
	mkdir -p build
	$(CC) $(CFLAGS) -o $@ $<

build/redelmeier: src/redelmeier.c
	mkdir -p build
	$(CC) $(CFLAGS) -o $@ $<

build/diag180: src/diag180.c FORCE180
	mkdir -p build
	$(CC) $(CFLAGS) $(D180FLAGS) -o $@ $<

test: build/redelmeier build/diag180
	sh tests/run_tests.sh
	sh tests/test_mirror.sh
	sh tests/test_diag180.sh

# diag180 is a CWEB program: src/diag180.w is the source.
# "make tangle" regenerates src/diag180.c (which is kept in the repository,
# so that ctangle is not needed for compiling);
# "make doc" produces doc/diag180.pdf (needs cweave, pdflatex and bibtex with
# the style plainurl; the bibliography is src/diag180.bib; the LaTeX
# class cweb-hy is in cweb/).
tangle:
	cd src && ctangle diag180.w - diag180.c

doc:
	mkdir -p doc
	cd doc && cweave ../src/diag180.w - diag180.tex
	cd doc && TEXINPUTS=../cweb: pdflatex -interaction=nonstopmode -halt-on-error diag180.tex >/dev/null
	cd doc && BIBINPUTS=../src: bibtex diag180 >/dev/null
	cd doc && TEXINPUTS=../cweb: pdflatex -interaction=nonstopmode -halt-on-error diag180.tex >/dev/null
	cd doc && TEXINPUTS=../cweb: pdflatex -interaction=nonstopmode -halt-on-error diag180.tex >/dev/null
	cd doc && TEXINPUTS=../cweb: pdflatex -interaction=nonstopmode -halt-on-error diag180.tex | grep -E "^Output|Warning.*undefined" ; true
	cd doc && rm -f diag180.aux diag180.log diag180.out diag180.toc diag180.idx diag180.scn diag180.tex diag180.bbl diag180.blg

# always rebuild diag180, so switching CNT128 on or off takes effect
FORCE180:

.PHONY: test FORCE180 tangle doc
