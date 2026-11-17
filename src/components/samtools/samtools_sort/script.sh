#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

[[ "$par_natural_sort" == "false" ]] && unset par_natural_sort
[[ "$par_ascii_sort" == "false" ]] && unset par_ascii_sort

cmd_args=(
  samtools sort
  ${meta_cpus:+-@ "$meta_cpus"}
  -T "$meta_temp_dir/samtools_sort_tmp"
  ${par_natural_sort:+-n}
  ${par_ascii_sort:+-N}
  ${par_compression:+-l "$par_compression"}
  -o "$par_output"
  "$par_input"
)

"${cmd_args[@]}"
