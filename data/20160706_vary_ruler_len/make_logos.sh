#!/usr/bin/env bash
set -euo pipefail

for i in {5..7}; do
    ../../tools/seq_logo -w300 mh_${i}_traj_??.tsv -o "mh_${i}_logo.pdf"
done
