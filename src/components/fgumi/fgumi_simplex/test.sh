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

grouped_bam="$meta_temp_dir/grouped.bam"
fgumi group --no-umi --strategy identity -i "$sorted_bam" -o "$grouped_bam"

log "Calling consensus with --min_reads 1"
out_bam="$meta_temp_dir/consensus.bam"
stats="$meta_temp_dir/stats.txt"
"$meta_executable" \
  --input "$grouped_bam" \
  --output "$out_bam" \
  --min_reads 1 \
  --read_name_prefix H0O_001 \
  --stats "$stats"

check_file_exists "$out_bam" "consensus BAM"
check_file_contains "$stats" "raw_reads_used" "simplex stats"

names=$(fgumi fastq -n -i "$out_bam" | awk 'NR % 4 == 1 { sub(/^@/, ""); print $1 }' | sort | uniq -c | awk '{ printf "%s:%s ", $2, $1 }')
if [[ "$names" != "H0O_001:0:2 H0O_001:1:2 H0O_001:2:2 H0O_001:3:2 H0O_001:4:2 " ]]; then
  log_error "unexpected consensus names (name:count): $names"
  exit 1
fi
log "✓ Five consensus pairs named H0O_001:0 to H0O_001:4"

log "Calling consensus with --min_reads 2"
out2_bam="$meta_temp_dir/consensus2.bam"
"$meta_executable" --input "$grouped_bam" --output "$out2_bam" --min_reads 2
n_pairs=$(fgumi fastq -n -i "$out2_bam" | awk 'NR % 4 == 1' | sort -u | wc -l)
[[ "$n_pairs" -eq 3 ]] || { log_error "expected 3 consensus pairs, got $n_pairs"; exit 1; }
log "✓ Three consensus pairs with --min_reads 2"

print_test_summary "fgumi_simplex tests"
