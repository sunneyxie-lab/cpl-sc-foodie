#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

unset_if_false=(
  par_remove_secondary
  par_disable_proper_pair_check
  par_add_template_cigar
  par_add_mate_score
  par_uncompressed
  par_no_PG
)
for par in "${unset_if_false[@]}"; do
  test_val="${!par}"
  [[ "$test_val" == "false" ]] && unset "$par"
done

cmd_args=(
  samtools fixmate
  ${meta_cpus:+-@ "$meta_cpus"}
  ${par_remove_secondary:+-r}
  ${par_disable_proper_pair_check:+-p}
  ${par_add_template_cigar:+-c}
  ${par_add_mate_score:+-m}
  ${par_uncompressed:+-u}
  ${par_sanitize:+-z "$par_sanitize"}
  ${par_no_PG:+--no-PG}
  ${par_output_fmt:+-O "$par_output_fmt"}
  "$par_input"
  "$par_output"
)

"${cmd_args[@]}"
