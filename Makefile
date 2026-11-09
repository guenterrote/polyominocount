CC = gcc
CFLAGS = -O3 -march=native -Wall

all: build/redelmeier build/ustats

build/ustats: src/ustats.c
	mkdir -p build
	$(CC) $(CFLAGS) -o $@ $<

build/redelmeier: src/redelmeier.c
	mkdir -p build
	$(CC) $(CFLAGS) -o $@ $<

test: build/redelmeier
	sh tests/run_tests.sh
	sh tests/test_mirror.sh

.PHONY: test
