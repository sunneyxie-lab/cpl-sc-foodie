#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

[[ "$par_non_directional" == "false" ]] && unset par_non_directional
[[ "$par_pbat" == "false" ]] && unset par_pbat

out="$meta_temp_dir/bismark_out"
tmp="$meta_temp_dir/bismark_tmp"
mkdir -p "$out"

# Each multicore chunk runs two bowtie2 instances in directional mode.
if [[ -n "$meta_cpus" && "$meta_cpus" -ge 2 ]]; then
  multicore=$(( meta_cpus / 2 ))
fi

cmd_args=(
  bismark
  --genome "$par_genome"
  -1 "$par_mates1"
  -2 "$par_mates2"
  -o "$out"
  -B bismark
  --temp_dir "$tmp"
  ${par_maxins:+-X "$par_maxins"}
  ${par_minins:+-I "$par_minins"}
  ${par_non_directional:+--non_directional}
  ${par_pbat:+--pbat}
  ${multicore:+--multicore "$multicore"}
)

"${cmd_args[@]}"

mv "$out/bismark_pe.bam" "$par_output"
mv "$out/bismark_PE_report.txt" "$par_report"
