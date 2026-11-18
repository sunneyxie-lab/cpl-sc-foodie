#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

IFS=';' read -r -a inputs <<< "$par_input"
out=$(mktemp "$meta_temp_dir/seqkit.XXXXXXXX")
cmd_args=(seqkit stats --all --tabular --threads 1 "${inputs[@]}")
"${cmd_args[@]}" > "$out"
mv "$out" "$par_out_file"
