#!/usr/bin/env zsh

git rev-parse HEAD > git_commit

for u in {5..7}; do
    for i in {01..05}; do
        ../../bin/mh                            \
            --num-moves 20000                   \
            --ruler-len $u                      \
            --temperature "300 5=>0"            \
            --random-seed $i                    \
            --output "mh_${u}_traj_${i}.tsv"    &
    done
done
wait
time
