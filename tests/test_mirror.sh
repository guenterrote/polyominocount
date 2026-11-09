#!/bin/sh
# Compare the transfer-matrix mirror counts with the baseline's mirror_v column.
set -e
N=${1:-14}
python3 proto/mirror_tm.py $N 2>/dev/null > /tmp/tm_$$.txt
./build/redelmeier -s $N 2>/dev/null | awk '/^# n fixed/{f=1;next} /^#/{f=0} f{print $1, $5}' > /tmp/rd_$$.txt
if diff /tmp/tm_$$.txt /tmp/rd_$$.txt; then echo "OK: mirror counts agree for n <= $N"; else echo MISMATCH; exit 1; fi
rm -f /tmp/tm_$$.txt /tmp/rd_$$.txt
