#!/bin/bash

## VIASH START
## VIASH END

source "$meta_resources_dir/test_helpers.sh"
setup_test_env

input_sam="$meta_temp_dir/input.sam"
seq=$(printf 'ACGTTGCAAC%.0s' 1 2 3 4 5)
qual=$(printf 'IIIIIIIIII%.0s' 1 2 3 4 5)

# Templates per (position, barcode): two cells share position 100, which must stay separate families.
printf '@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:10000\n' > "$input_sam"
i=0
for rec in "100 AAAA 3" "100 CCCC 2" "1000 AAAA 1" "2000 CCCC 2" "3000 GGGG 1"; do
  set -- $rec
  mpos=$(($1 + 150))
  for _ in $(seq 1 "$3"); do
    i=$((i + 1))
    printf 'q%s\t99\tchr1\t%s\t60\t50M\t=\t%s\t200\t%s\t%s\tMC:Z:50M\tCB:Z:%s\n' "$i" "$1" "$mpos" "$seq" "$qual" "$2" >> "$input_sam"
    printf 'q%s\t147\tchr1\t%s\t60\t50M\t=\t%s\t-200\t%s\t%s\tMC:Z:50M\tCB:Z:%s\n' "$i" "$mpos" "$1" "$seq" "$qual" "$2" >> "$input_sam"
  done
done

sorted_bam="$meta_temp_dir/sorted.bam"
fgumi sort --order template-coordinate -i "$input_sam" -o "$sorted_bam"

log "Grouping by position and cell barcode"
out_bam="$meta_temp_dir/grouped.bam"
fam="$meta_temp_dir/fam.txt"
grp="$meta_temp_dir/grp.txt"
"$meta_executable" \
  --input "$sorted_bam" \
  --output "$out_bam" \
  --strategy identity \
  --no_umi \
  --family_size_histogram "$fam" \
  --grouping_metrics "$grp"

check_file_exists "$out_bam" "grouped BAM"
check_file_not_empty "$out_bam" "grouped BAM"
check_file_exists "$fam" "family size histogram"
check_file_exists "$grp" "grouping metrics"

hist=$(awk 'NR > 1 { printf "%s:%s ", $1, $2 }' "$fam")
if [[ "$hist" != "1:2 2:2 3:1 " ]]; then
  log_error "unexpected family size histogram (size:count): $hist"
  exit 1
fi
log "✓ Families split by cell barcode: $hist"

accepted=$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "accepted_sam_records") c = i; next } { print $c }' "$grp")
[[ "$accepted" -eq 18 ]] || { log_error "expected 18 accepted records, got $accepted"; exit 1; }
log "✓ All 18 records accepted"

print_test_summary "fgumi_group tests"
