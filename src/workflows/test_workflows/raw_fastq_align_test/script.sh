#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

cmd=(
  python3 "$meta_resources_dir/check.py"
  --id "$par_id"
  --read1 "$par_read1" --read2 "$par_read2"
  --output_bam "$par_output_bam"
  --output_bismark_report "$par_output_bismark_report"
  --output_raw_stats "$par_output_raw_stats"
  --output_extracted_trimmed_stats "$par_output_extracted_trimmed_stats"
  --output_raw_fastqc_html "$par_output_raw_fastqc_html"
  --output_raw_fastqc_zip "$par_output_raw_fastqc_zip"
  --output_extracted_fastqc_html "$par_output_extracted_fastqc_html"
  --output_extracted_fastqc_zip "$par_output_extracted_fastqc_zip"
  --output_trimmed_fastqc_html "$par_output_trimmed_fastqc_html"
  --output_trimmed_fastqc_zip "$par_output_trimmed_fastqc_zip"
)
"${cmd[@]}"
