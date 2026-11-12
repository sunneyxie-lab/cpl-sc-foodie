#!/bin/bash
# All paths come from the login driver: Slurm spools this script elsewhere.
set -euo pipefail

if [[ $# -ne 6 ]]; then
  echo 'usage: common_extract.sh <reference> <chr14.fa> <source-list-1> <list-1> <source-list-2> <list-2>' >&2
  exit 1
fi
for path in "$@"; do
  [[ "$path" == /* ]] || { echo "error: expected absolute path: $path" >&2; exit 1; }
done

# Consume the full stream; an early awk exit would SIGPIPE zcat under pipefail.
zcat -f -- "$1" | awk '
  /^>/ { split(substr($0, 2), words, /[ \t]/); keep = (words[1] == "chr14") }
  keep
' > "$2"
[[ -s "$2" ]] || { echo "error: chr14 is absent from $1" >&2; exit 1; }
cp -- "$3" "$4"
cp -- "$5" "$6"
[[ -s "$4" && -s "$6" ]] || { echo 'error: copied cell list is empty' >&2; exit 1; }
