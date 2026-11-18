#!/bin/bash

## VIASH START
## VIASH END

source "$meta_resources_dir/test_helpers.sh"
setup_test_env

input_sam="$meta_temp_dir/input.sam"
seq=$(printf 'ACGTTGCAAC%.0s' 1 2 3 4 5)
qual=$(printf 'IIIIIIIIII%.0s' 1 2 3 4 5)

printf '@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:10000\n' > "$input_sam"
i=0
for rec in "3000 GGGG" "100 AAAA" "2000 CCCC" "1000 AAAA" "500 CCCC"; do
  set -- $rec
  i=$((i + 1))
  mpos=$(($1 + 150))
  printf 'q%s\t99\tchr1\t%s\t60\t50M\t=\t%s\t200\t%s\t%s\tMC:Z:50M\tCB:Z:%s\n' "$i" "$1" "$mpos" "$seq" "$qual" "$2" >> "$input_sam"
  printf 'q%s\t147\tchr1\t%s\t60\t50M\t=\t%s\t-200\t%s\t%s\tMC:Z:50M\tCB:Z:%s\n' "$i" "$mpos" "$1" "$seq" "$qual" "$2" >> "$input_sam"
done

log "Sorting into template-coordinate order"
out_bam="$meta_temp_dir/out.bam"
"$meta_executable" --input "$input_sam" --output "$out_bam" --order template-coordinate

check_file_exists "$out_bam" "sorted BAM"
magic=$(head -c2 "$out_bam" | od -An -tx1)
[[ "$magic" == *"1f 8b"* ]] || { log_error "output is not gzip/BGZF: $magic"; exit 1; }
fgumi sort --verify --order template-coordinate -i "$out_bam" || { log_error "output is not template-coordinate sorted"; exit 1; }
log "✓ Output is a template-coordinate sorted BAM"

n_reads=$(fgumi fastq -i "$out_bam" | awk 'NR % 4 == 1' | wc -l)
[[ "$n_reads" -eq 10 ]] || { log_error "expected 10 reads, got $n_reads"; exit 1; }
log "✓ All 10 reads kept"

print_test_summary "fgumi_sort tests"
