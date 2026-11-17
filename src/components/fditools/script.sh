#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

out=$(mktemp -d "$meta_temp_dir/fditools.XXXXXXXX")
cpus=${meta_cpus:-1}

fdi --version
cmd_args=(
  fdi cell-barcode
  -a "$par_read1" -b "$par_read2"
  -c "$out/read1.fastq.gz" -d "$out/read2.fastq.gz"
  -m 11 -s 19 -t "$cpus" -w "$par_bc_whitelist"
  -T 1 -Z
)
printf 'Running:' >&2
printf ' %q' "${cmd_args[@]}" >&2
printf '\n' >&2
"${cmd_args[@]}"

mv "$out/read1.fastq.gz" "$par_output_read1"
mv "$out/read2.fastq.gz" "$par_output_read2"
