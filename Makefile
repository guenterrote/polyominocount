CC = gcc
CFLAGS = -O3 -march=native -Wall

all: build/redelmeier build/ustats build/diag180

build/ustats: src/ustats.c
	mkdir -p build
	$(CC) $(CFLAGS) -o $@ $<

build/redelmeier: src/redelmeier.c
	mkdir -p build
	$(CC) $(CFLAGS) -o $@ $<

build/diag180: src/diag180.c
	mkdir -p build
	$(CC) $(CFLAGS) -o $@ $<

test: build/redelmeier build/diag180
	sh tests/run_tests.sh
	sh tests/test_mirror.sh
	sh tests/test_diag180.sh

.PHONY: test
