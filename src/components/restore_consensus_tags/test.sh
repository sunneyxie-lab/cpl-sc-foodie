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
  "H0O_001:0	77	*	0	0	*	*	0	0	$seq	$qual	MI:Z:0	CB:Z:AAAA	cD:i:3	cM:i:3" \
  "H0O_001:0	141	*	0	0	*	*	0	0	$seq	$qual	MI:Z:0	CB:Z:AAAA	cD:i:2	cM:i:2" \
  "H0O_001:1	77	*	0	0	*	*	0	0	$seq	$qual	MI:Z:1	cD:i:1	CB:Z:CCCC" \
  "H0O_001:1	141	*	0	0	*	*	0	0	$seq	$qual	MI:Z:1	cD:i:1	CB:Z:CCCC" \
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

log "Restoring consensus tags"
out_bam="$meta_temp_dir/out.bam"
"$meta_executable" --input "$aligned_bam" --consensus "$consensus_bam" --output "$out_bam"

check_file_exists "$out_bam" "tagged BAM"
samtools view "$out_bam" | awk -F'\t' '
  {
    expected = ($1 ~ /^H0O_001:0_/) ? "AAAA" : "CCCC"
    split($1, parts, "_CB:Z:")
    if (parts[2] != expected || $(NF - 1) != "CB:Z:" expected) exit 1
  }
' || { log_error "read name suffix or CB tag is wrong"; exit 1; }
log "✓ Every record is renamed <name>_CB:Z:<bc> and tagged CB:Z:<bc>"

actual_cd=$(samtools view "$out_bam" | awk -F'\t' '{ print $NF }' | tr '\n' ' ')
[[ "$actual_cd" == "cD:i:3 cD:i:2 cD:i:1 cD:i:1 " ]] || { log_error "unexpected cD tags: $actual_cd"; exit 1; }
log "✓ Each record gets the cD of its consensus read (same name, same R1/R2)"

samtools view "$out_bam" | awk -F'\t' -v OFS='\t' '{ sub(/_CB:Z:[ACGT]+$/, "", $1); NF -= 2; print }' > "$meta_temp_dir/stripped.sam"
samtools view "$aligned_bam" > "$meta_temp_dir/expected.sam"
cmp "$meta_temp_dir/stripped.sam" "$meta_temp_dir/expected.sam" || { log_error "records changed beyond name, CB and cD tags"; exit 1; }
log "✓ Records are otherwise unchanged"

log "Checking that a read without a consensus mate fails the run"
missing_bam="$meta_temp_dir/missing.bam"
{ cat "$aligned_sam"; printf 'H0O_001:9\t99\tchr1\t900\t42\t10M\t=\t950\t60\t%s\t%s\n' "$seq" "$qual"; } |
  samtools view -b -o "$missing_bam" -
if "$meta_executable" --input "$missing_bam" --consensus "$consensus_bam" --output "$meta_temp_dir/fail.bam"; then
  log_error "expected a non-zero exit for a read without a consensus mate"
  exit 1
fi
log "✓ A read without a consensus mate fails the run"

log "Checking that a consensus read without cD fails the run"
no_cd_bam="$meta_temp_dir/no_cd.bam"
awk -F'\t' -v OFS='\t' '$1 == "H0O_001:1" && $2 == 141 { $13 = "RG:Z:x" } { print }' "$consensus_sam" |
  samtools view -b -o "$no_cd_bam" -
if "$meta_executable" --input "$aligned_bam" --consensus "$no_cd_bam" --output "$meta_temp_dir/fail_cd.bam"; then
  log_error "expected a non-zero exit for a consensus read without cD"
  exit 1
fi
log "✓ A consensus read without cD fails the run"

print_test_summary "restore_consensus_tags tests"
