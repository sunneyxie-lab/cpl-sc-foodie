#!/bin/bash

## VIASH START
## VIASH END

source "$meta_resources_dir/test_helpers.sh"
setup_test_env

seq=ACGTACGTAC
qual=IIIIIIIIII

consensus_sam="$meta_temp_dir/consensus.sam"
consensus_bam="$meta_temp_dir/consensus.bam"
printf '@HD\tVN:1.6\tSO:unsorted\n' > "$consensus_sam"
printf '%s\n' \
  "H0O_001:0	77	*	0	0	*	*	0	0	$seq	$qual	MI:Z:0	CB:Z:AAAA" \
  "H0O_001:0	141	*	0	0	*	*	0	0	$seq	$qual	MI:Z:0	CB:Z:AAAA" \
  "H0O_001:1	77	*	0	0	*	*	0	0	$seq	$qual	MI:Z:1	CB:Z:CCCC" \
  "H0O_001:1	141	*	0	0	*	*	0	0	$seq	$qual	MI:Z:1	CB:Z:CCCC" \
  >> "$consensus_sam"
samtools view -b -o "$consensus_bam" "$consensus_sam"

aligned_sam="$meta_temp_dir/aligned.sam"
aligned_bam="$meta_temp_dir/aligned.bam"
printf '@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:10000\n' > "$aligned_sam"
printf '%s\n' \
  "H0O_001:0	99	chr1	100	42	10M	=	300	210	$seq	$qual	NM:i:0	XR:Z:CT	XG:Z:CT" \
  "H0O_001:0	147	chr1	300	42	10M	=	100	-210	$seq	$qual	NM:i:0	XR:Z:GA	XG:Z:CT" \
  "H0O_001:1	83	chr1	700	42	10M	=	500	-210	$seq	$qual	NM:i:1	XR:Z:CT	XG:Z:GA" \
  "H0O_001:1	163	chr1	500	42	10M	=	700	210	$seq	$qual	NM:i:0	XR:Z:GA	XG:Z:GA" \
  >> "$aligned_sam"
samtools view -b -o "$aligned_bam" "$aligned_sam"

log "Restoring cell barcodes"
out_bam="$meta_temp_dir/out.bam"
"$meta_executable" --input "$aligned_bam" --consensus "$consensus_bam" --output "$out_bam"

check_file_exists "$out_bam" "tagged BAM"
samtools view "$out_bam" | awk -F'\t' '
  {
    expected = ($1 ~ /^H0O_001:0_/) ? "AAAA" : "CCCC"
    split($1, parts, "_CB:Z:")
    if (parts[2] != expected || $NF != "CB:Z:" expected) exit 1
  }
' || { log_error "read name suffix or CB tag is wrong"; exit 1; }
log "✓ Every record is renamed <name>_CB:Z:<bc> and tagged CB:Z:<bc>"

samtools view "$out_bam" | awk -F'\t' -v OFS='\t' '{ sub(/_CB:Z:[ACGT]+$/, "", $1); NF--; print }' > "$meta_temp_dir/stripped.sam"
samtools view "$aligned_bam" > "$meta_temp_dir/expected.sam"
cmp "$meta_temp_dir/stripped.sam" "$meta_temp_dir/expected.sam" || { log_error "records changed beyond name and CB tag"; exit 1; }
log "✓ Records are otherwise unchanged"

log "Checking that a read without a barcode fails the run"
missing_bam="$meta_temp_dir/missing.bam"
{ cat "$aligned_sam"; printf 'H0O_001:9\t99\tchr1\t900\t42\t10M\t=\t950\t60\t%s\t%s\n' "$seq" "$qual"; } |
  samtools view -b -o "$missing_bam" -
if "$meta_executable" --input "$missing_bam" --consensus "$consensus_bam" --output "$meta_temp_dir/fail.bam"; then
  log_error "expected a non-zero exit for a read without a barcode"
  exit 1
fi
log "✓ A read without a barcode fails the run"

print_test_summary "restore_cell_barcode tests"
