#!/bin/bash

## VIASH START
## VIASH END

source "$meta_resources_dir/test_helpers.sh"
setup_test_env

seq=ACGTACGTAC
qual=IIIIIIIIII

input_sam="$meta_temp_dir/input.sam"
input_bam="$meta_temp_dir/input.bam"
printf '@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:10000\n' > "$input_sam"
printf '%s\n' \
  "LH1:8:2466:17478_1:N:0:GAACTGAG+TACTGTTC_CB:Z:TATTGTGG	99	chr1	100	60	10M	=	300	210	$seq	$qual	NM:i:0	XG:Z:CT" \
  "LH1:8:2466:17478_1:N:0:GAACTGAG+TACTGTTC_CB:Z:TATTGTGG	147	chr1	300	60	10M	=	100	-210	$seq	$qual	NM:i:0	XG:Z:CT" \
  "odd_CB:Z:x_name_CB:Z:CCCCAAAA	83	chr1	700	60	10M	=	500	-210	$seq	$qual" \
  "odd_CB:Z:x_name_CB:Z:CCCCAAAA	163	chr1	500	60	10M	=	700	210	$seq	$qual" \
  >> "$input_sam"
samtools view -b -o "$input_bam" "$input_sam"

log "Tagging cell barcodes"
out_bam="$meta_temp_dir/out.bam"
"$meta_executable" --input "$input_bam" --output "$out_bam"

check_file_exists "$out_bam" "tagged BAM"
tags=$(samtools view "$out_bam" | awk -F'\t' '{ printf "%s ", $NF }')
if [[ "$tags" != "CB:Z:TATTGTGG CB:Z:TATTGTGG CB:Z:CCCCAAAA CB:Z:CCCCAAAA " ]]; then
  log_error "unexpected CB tags (barcode taken after the last _CB:Z:): $tags"
  exit 1
fi
log "✓ Every record gets the barcode after its last _CB:Z:"

samtools view "$out_bam" | awk -F'\t' -v OFS='\t' '{ NF--; print }' > "$meta_temp_dir/stripped.sam"
samtools view "$input_bam" > "$meta_temp_dir/expected.sam"
cmp "$meta_temp_dir/stripped.sam" "$meta_temp_dir/expected.sam" || { log_error "records changed beyond the added tag"; exit 1; }
log "✓ Records are otherwise unchanged"

log "Checking that a read without a barcode suffix fails the run"
bad_bam="$meta_temp_dir/bad.bam"
{ cat "$input_sam"; printf 'plain_read\t99\tchr1\t900\t60\t10M\t=\t950\t60\t%s\t%s\n' "$seq" "$qual"; } |
  samtools view -b -o "$bad_bam" -
if "$meta_executable" --input "$bad_bam" --output "$meta_temp_dir/fail.bam"; then
  log_error "expected a non-zero exit for a read without a barcode suffix"
  exit 1
fi
log "✓ A read without a barcode suffix fails the run"

print_test_summary "tag_cell_barcode tests"
