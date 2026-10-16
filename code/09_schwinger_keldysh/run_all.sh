#!/bin/sh
# Run every check of this folder and store the outputs in expected_output/.
# Usage:  cd code/09_schwinger_keldysh && sh run_all.sh      (total runtime about 2 minutes)
set -e
cd "$(dirname "$0")"
mkdir -p expected_output
for s in tables234 numcheck table5 numcheck_k analysis2 analysis3 extract; do
  echo "== $s"; python3 $s.py > expected_output/$s.txt 2>&1
done
for n in 2 4 6 10; do
  echo "== kexp $n"; python3 kexp.py $n > expected_output/kexp_N$n.txt 2>&1
done
echo "done"
