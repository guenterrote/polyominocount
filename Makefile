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

# always rebuild diag180, so switching CNT128 on or off takes effect
FORCE180:

.PHONY: test FORCE180
