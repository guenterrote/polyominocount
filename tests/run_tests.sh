#!/bin/sh
# Compare the C counter (fixed and free columns) with the brute-force Python reference.
set -e
N=${1:-10}
python3 tests/brute.py $N > /tmp/brute_$$.txt
./build/redelmeier -s $N 2>/dev/null | awk '/^# n fixed/{f=1;next} /^#/{f=0} f{print $1, $2, $NF}' > /tmp/c_$$.txt
if diff /tmp/brute_$$.txt /tmp/c_$$.txt; then echo "OK: fixed and free counts agree for n <= $N"; else echo "MISMATCH"; exit 1; fi
rm -f /tmp/brute_$$.txt /tmp/c_$$.txt
