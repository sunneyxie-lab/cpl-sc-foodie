#!/bin/bash

## VIASH START
## VIASH END

source "$meta_resources_dir/test_helpers.sh"
setup_test_env

seq=ACGTACGTAC
qual=IIIIIIIIII

input_sam="$meta_temp_dir/input.sam"
input_bam="$meta_temp_dir/input.bam"
printf '@HD\tVN:1.0\tSO:unsorted\tGO:query\tSS:unsorted:template-coordinate\n@SQ\tSN:chr1\tLN:10000\n@PG\tID:fgumi\tPN:fgumi\n' > "$input_sam"
i=0
for cb in AAAA CCCC GGGG AAAA GGGG CCCC AAAA GGGG; do
  mi=$((i / 2))
  printf 'r%d\t99\tchr1\t%d\t60\t10M\t=\t%d\t210\t%s\t%s\tCB:Z:%s\tMI:Z:%d\n' "$i" $((100 + i)) $((300 + i)) "$seq" "$qual" "$cb" "$mi" >> "$input_sam"
  printf 'r%d\t147\tchr1\t%d\t60\t10M\t=\t%d\t-210\t%s\t%s\tCB:Z:%s\tMI:Z:%d\n' "$i" $((300 + i)) $((100 + i)) "$seq" "$qual" "$cb" "$mi" >> "$input_sam"
  i=$((i + 1))
done
samtools view -b -o "$input_bam" "$input_sam"

cell_list="$meta_temp_dir/cell_list.tsv"
printf 'AAAA\t5\nGGGG\t3\n' > "$cell_list"

log "Filtering to the listed cells"
out_bam="$meta_temp_dir/out.bam"
"$meta_executable" --input "$input_bam" --cell_list "$cell_list" --output "$out_bam"

samtools quickcheck "$out_bam" || { log_error "output fails samtools quickcheck"; exit 1; }
log "✓ Output is a valid BAM"

samtools view "$out_bam" > "$meta_temp_dir/got.sam"
samtools view "$input_bam" |
  awk -F'\t' '{ for (i = 12; i <= NF; i++) if ($i == "CB:Z:AAAA" || $i == "CB:Z:GGGG") { print; break } }' \
  > "$meta_temp_dir/expected.sam"
[[ -s "$meta_temp_dir/expected.sam" ]] || { log_error "test input has no records of listed cells"; exit 1; }
cmp "$meta_temp_dir/got.sam" "$meta_temp_dir/expected.sam" || { log_error "records differ from the listed cells' input records"; exit 1; }
log "✓ Output has exactly the listed cells' records, in input order"

samtools view --no-PG -H "$input_bam" > "$meta_temp_dir/in.header"
samtools view --no-PG -H "$out_bam" > "$meta_temp_dir/out.header"
[[ "$(grep '^@HD' "$meta_temp_dir/out.header")" == "$(grep '^@HD' "$meta_temp_dir/in.header")" ]] ||
  { log_error "@HD line changed: $(grep '^@HD' "$meta_temp_dir/out.header")"; exit 1; }
log "✓ @HD line is unchanged"
cmp <(grep -v '^@PG' "$meta_temp_dir/out.header") <(grep -v '^@PG' "$meta_temp_dir/in.header") ||
  { log_error "header changed beyond @PG lines"; exit 1; }
n_pg_in=$(grep -c '^@PG' "$meta_temp_dir/in.header")
n_pg_out=$(grep -c '^@PG' "$meta_temp_dir/out.header")
[[ "$n_pg_out" -eq $((n_pg_in + 1)) ]] || { log_error "expected one added @PG line, got $n_pg_in -> $n_pg_out"; exit 1; }
log "✓ Header differs only by one added @PG line"

log "Checking that a cell list matching no record fails the run"
printf 'TTTT\t1\n' > "$meta_temp_dir/other_cells.tsv"
if "$meta_executable" --input "$input_bam" --cell_list "$meta_temp_dir/other_cells.tsv" --output "$meta_temp_dir/none.bam"; then
  log_error "expected a non-zero exit for a cell list matching no record"
  exit 1
fi
log "✓ A cell list matching no record fails the run"

print_test_summary "filter_cells tests"
