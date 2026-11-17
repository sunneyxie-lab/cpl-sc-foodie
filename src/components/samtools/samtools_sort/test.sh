#!/bin/bash

## VIASH START
## VIASH END

source "$meta_resources_dir/test_helpers.sh"
setup_test_env

input_sam="$meta_temp_dir/input.sam"
input_bam="$meta_temp_dir/input.bam"

printf '@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:10000\n' > "$input_sam"
for rec in "read10 500" "read2 100" "read1 900" "read3 300"; do
  set -- $rec
  printf '%s\t0\tchr1\t%s\t60\t4M\t*\t0\t0\tACGT\tIIII\n' "$1" "$2" >> "$input_sam"
done
samtools view -b -o "$input_bam" "$input_sam"

log "Sorting by read name (natural)"
name_bam="$meta_temp_dir/name.bam"
"$meta_executable" --input "$input_bam" --output "$name_bam" --natural_sort

check_file_exists "$name_bam" "name-sorted BAM"
samtools view -H "$name_bam" | grep -q '^@HD.*SO:queryname' || { log_error "name-sorted BAM header lacks SO:queryname"; exit 1; }
names=$(samtools view "$name_bam" | cut -f1 | tr '\n' ' ')
if [[ "$names" != "read1 read2 read3 read10 " ]]; then
  log_error "Unexpected name order: $names"
  exit 1
fi
log "✓ Read names come out in natural order"
check_file_not_exists "$name_bam.txt" "SAM text dump"

log "Sorting by coordinate"
coord_bam="$meta_temp_dir/coord.bam"
"$meta_executable" --input "$input_bam" --output "$coord_bam" --compression 6

check_file_exists "$coord_bam" "coordinate-sorted BAM"
samtools view -H "$coord_bam" | grep -q '^@HD.*SO:coordinate' || { log_error "coordinate-sorted BAM header lacks SO:coordinate"; exit 1; }
samtools view "$coord_bam" | awk '$4 < prev { exit 1 } { prev = $4 }' || { log_error "positions are not non-decreasing"; exit 1; }
[[ $(samtools view -c "$coord_bam") -eq 4 ]] || { log_error "record count changed"; exit 1; }
log "✓ Positions are non-decreasing"
check_file_not_exists "$coord_bam.txt" "SAM text dump"

print_test_summary "samtools_sort tests"
