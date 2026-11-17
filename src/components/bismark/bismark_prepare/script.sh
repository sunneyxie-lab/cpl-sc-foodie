#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

mkdir -p "$par_output"

# Copies, not symlinks: bismark align reads the FASTA from this directory, and
# rclone skips symlinks when the directory is uploaded as a test resource.
IFS=';' read -ra inputs <<< "$par_input"
for fasta in "${inputs[@]}"; do
  cp "$fasta" "$par_output/"
done

# Two indexers run at once, each with --parallel threads.
if [[ -n "$meta_cpus" && "$meta_cpus" -ge 4 ]]; then
  parallel=$(( meta_cpus / 2 ))
fi

cmd_args=(
  bismark prepare
  ${parallel:+--parallel "$parallel"}
  "$par_output"
)

"${cmd_args[@]}"
