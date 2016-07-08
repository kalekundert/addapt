#!/usr/bin/env zsh

git rev-parse HEAD > git_commit

for i in {01..05}; do
    ../../bin/mh -n 20000 -r $i -T "300 5=>0" -o "mh_$i.tsv" &
done
wait
time
