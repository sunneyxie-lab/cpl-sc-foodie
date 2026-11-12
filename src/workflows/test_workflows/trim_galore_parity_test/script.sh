#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

work_dir=$(mktemp -d "$meta_temp_dir/trim_galore_parity.XXXXXX")
cmd=(
  python3 "$meta_resources_dir/compare.py"
  --id "$par_id"
  --extracted_read1 "$par_extracted_read1"
  --extracted_read2 "$par_extracted_read2"
  --rust_read1 "$par_rust_read1" --rust_read2 "$par_rust_read2"
  --work_dir "$work_dir"
  --output "$work_dir/evidence.json"
)
"${cmd[@]}"
mv "$work_dir/evidence.json" "$par_output"
