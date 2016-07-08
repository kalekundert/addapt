#!/usr/bin/env zsh

git rev-parse HEAD > git_commit

for i in {01..04}; do
    ../../bin/mh                            \
        --num-moves 40000                   \
        --temperature "500 5=>0"            \
        --output "run_$i.tsv"               &
done
wait
time
