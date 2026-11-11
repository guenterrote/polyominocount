#!/bin/sh
# Compare fix(r180) from the diagonal transfer matrix (build/diag180) with the
# r180 column of the Redelmeier baseline (build/redelmeier -s).
set -e
N=${1:-16}
./build/redelmeier -s $N 2>/dev/null | awk '/^# n fixed/{f=1;next} /^#/{f=0} f{print $1, $4}' > /tmp/r_$$.txt
./build/diag180 -o /tmp/res_$$.py $N 2>/dev/null | grep -v "^#" | awk '{print $1, $2}' > /tmp/d_$$.txt
if diff /tmp/r_$$.txt /tmp/d_$$.txt; then echo "OK: diag180 agrees with Redelmeier r180 for n <= $N"; else echo "MISMATCH"; exit 1; fi
./build/diag180 -o /tmp/res_$$.py -m 1000003 $N 2>/dev/null | grep -v "^#" | awk 'NR==FNR{a[$1]=$2%1000003;next} a[$1]!=$2{bad=1} END{exit bad}' /tmp/d_$$.txt - && echo "OK: modular counters agree"
# the production version (src/production.ch) must give the same table
./build/diag180 -v 0 -o /tmp/res_$$.py 24 2>/dev/null | grep -v "^#" > /tmp/d_$$.txt
./build/diag180p -v 0 -o /tmp/res_$$.py 24 2>/dev/null | grep -v "^#" > /tmp/p_$$.txt
if cmp -s /tmp/d_$$.txt /tmp/p_$$.txt && [ -s /tmp/p_$$.txt ]; then echo "OK: production version diag180p agrees with diag180 for n <= 24"; else echo "MISMATCH diag180p"; exit 1; fi
rm -f /tmp/p_$$.txt
rm -f /tmp/r_$$.txt /tmp/d_$$.txt /tmp/res_$$.py
