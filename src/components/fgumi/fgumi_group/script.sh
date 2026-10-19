#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

[[ "$par_no_umi" == "false" ]] && unset par_no_umi

cmd_args=(
  fgumi group
  -i "$par_input"
  -o "$par_output"
  ${par_family_size_histogram:+-f "$par_family_size_histogram"}
  ${par_grouping_metrics:+-g "$par_grouping_metrics"}
  -s "$par_strategy"
  ${par_no_umi:+--no-umi}
  ${par_min_map_q:+-m "$par_min_map_q"}
  ${meta_cpus:+--threads "$meta_cpus"}
)

"${cmd_args[@]}"
