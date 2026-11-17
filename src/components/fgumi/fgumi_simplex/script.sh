#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

cmd_args=(
  fgumi simplex
  -i "$par_input"
  -o "$par_output"
  ${par_stats:+-s "$par_stats"}
  -M "$par_min_reads"
  ${par_read_name_prefix:+-p "$par_read_name_prefix"}
  ${meta_cpus:+--threads "$meta_cpus"}
)

"${cmd_args[@]}"
