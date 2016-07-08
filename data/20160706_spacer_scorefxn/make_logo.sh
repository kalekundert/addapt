#!/usr/bin/env bash
set -euo pipefail

../../tools/seq_logo -w350 run_??.tsv -o "logo.pdf"
