#!/bin/bash

## VIASH START
## VIASH END

source "$meta_resources_dir/test_helpers.sh"
setup_test_env

input_sam="$meta_temp_dir/input.sam"
input_bam="$meta_temp_dir/input.bam"

printf '@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:10000\n' > "$input_sam"
printf '%s\n' \
  "p1	99	chr1	100	60	10M	=	300	210	ACGTACGTAC	IIIIIIIIII	CB:Z:AAAA" \
  "p1	147	chr1	300	60	10M	=	100	-210	ACGTACGTAC	IIIIIIIIII	CB:Z:AAAA" \
  "p2	99	chr1	500	60	10M	=	700	210	ACGTACGTAC	IIIIIIIIII	CB:Z:CCCC" \
  "p2	147	chr1	700	60	10M	=	500	-210	ACGTACGTAC	IIIIIIIIII	CB:Z:CCCC" \
  >> "$input_sam"
samtools view -b -o "$input_bam" "$input_sam"

log "Running fixmate with defaults"
out_bam="$meta_temp_dir/out.bam"
"$meta_executable" --input "$input_bam" --output "$out_bam"

check_file_exists "$out_bam" "fixmate BAM"
[[ $(samtools view -c "$out_bam") -eq 4 ]] || { log_error "expected 4 records"; exit 1; }
samtools view "$out_bam" | awk '!/\tMC:Z:/ { exit 1 }' || { log_error "a record lacks MC:Z:"; exit 1; }
samtools view "$out_bam" | awk '!/\tCB:Z:(AAAA|CCCC)/ { exit 1 }' || { log_error "a record lost its CB:Z: tag"; exit 1; }
log "✓ Every record has MC:Z: and keeps CB:Z:"

log "Running fixmate with --add_mate_score"
ms_bam="$meta_temp_dir/ms.bam"
"$meta_executable" --input "$input_bam" --output "$ms_bam" --add_mate_score
samtools view "$ms_bam" | awk '!/\tms:i:/ { exit 1 }' || { log_error "a record lacks ms:i:"; exit 1; }
log "✓ Every record has ms:i:"

print_test_summary "samtools_fixmate tests"
