#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# -D wants bare barcodes, one per line, not the two-column cell list.
cut -f1 "$par_cell_list" > "$meta_temp_dir/cell_barcodes.txt"
samtools view -b ${meta_cpus:+-@ "$meta_cpus"} -D "CB:$meta_temp_dir/cell_barcodes.txt" -o "$par_output" "$par_input"

if [[ $(samtools view -c ${meta_cpus:+-@ "$meta_cpus"} "$par_output") -eq 0 ]]; then
  echo "filter_cells: no record of $par_input has a CB in $par_cell_list" >&2
  exit 1
fi
