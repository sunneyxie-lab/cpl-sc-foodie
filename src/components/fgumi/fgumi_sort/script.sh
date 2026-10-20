#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

cmd_args=(
  fgumi sort
  -i "$par_input"
  -o "$par_output"
  ${par_order:+--order "$par_order"}
  ${par_key_types:+--key-types "$par_key_types"}
  ${meta_cpus:+--threads "$meta_cpus"}
  --tmp-dir "$meta_temp_dir"
)

"${cmd_args[@]}"
