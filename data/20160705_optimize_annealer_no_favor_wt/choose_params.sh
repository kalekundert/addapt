#!/usr/bin/env bash
set -euo pipefail

MH=../../bin/mh

git rev-parse HEAD > git_commit

for N in {200,300,400,500,600}; do
    $MH -n 2000 -T "$N 5=>0" -o "N_$N.tsv" &
done

for T in {1,2,3,4,5,6}; do
    $MH -n 2000 -T "$((300 * T / 5)) $T=>0" -o "T_hi_$T.tsv" &
done

wait


