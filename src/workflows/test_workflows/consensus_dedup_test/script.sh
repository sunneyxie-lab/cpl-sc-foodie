#!/bin/bash

## VIASH START
## VIASH END

source "$meta_resources_dir/test_helpers.sh"
setup_test_env

samtools quickcheck "$par_output_bam" || { log_error "BAM fails samtools quickcheck: $par_output_bam"; exit 1; }
n_records=$(samtools view -c "$par_output_bam")
[[ "$n_records" -gt 0 ]] || { log_error "BAM has no records"; exit 1; }
log "✓ BAM is valid and has $n_records records"

samtools view -H "$par_output_bam" | grep -q '^@HD.*SO:coordinate' || { log_error "BAM header lacks SO:coordinate"; exit 1; }
log "✓ BAM is coordinate-sorted"

check_file_not_empty "$par_output_bai" "BAM index"
bai_magic=$(head -c4 "$par_output_bai" | od -An -tx1 | tr -d ' ')
[[ "$bai_magic" == "42414901" ]] || { log_error "index does not start with BAI\\1: $bai_magic"; exit 1; }
log "✓ Index is a BAI file"

check_file_not_empty "$par_bismark_report" "Bismark report"
check_file_not_empty "$par_family_size_histogram" "family size histogram"
check_file_not_empty "$par_grouping_metrics" "grouping metrics"
check_file_not_empty "$par_simplex_stats" "simplex stats"
check_file_not_empty "$par_cell_list" "cell list"
if [[ -n "${par_input_cell_list:-}" ]]; then
  cmp "$par_cell_list" "$par_input_cell_list" || { log_error "published cell list differs from the passed-in one"; exit 1; }
  log "✓ Published cell list equals the passed-in list"
  if [[ -n "${par_cell_counts:-}${par_cell_qc:-}${par_knee_plot:-}" ]]; then
    log_error "cell-calling outputs were published although a cell list was passed in"
    exit 1
  fi
  log "✓ No cell-calling outputs with a passed-in list"
else
  check_file_not_empty "$par_cell_counts" "cell counts"
  check_file_not_empty "$par_cell_qc" "cell QC summary"
  check_file_not_empty "$par_knee_plot" "knee plot"
fi

out_barcodes="$meta_temp_dir/out_barcodes.txt"
# IDs may contain "_", so split at the last "_CB:Z:" and compare the ID as a literal string.
samtools view "$par_output_bam" | awk -F'\t' -v id="$par_id" -v out="$out_barcodes" '
  {
    name = $1; cut = 0; rest = name
    while ((k = index(rest, "_CB:Z:")) > 0) { cut += k; rest = substr(rest, k + 1) }
    if (cut == 0) { print "read name lacks _CB:Z: suffix: " name > "/dev/stderr"; exit 1 }
    prefix = substr(name, 1, cut - 1)
    suffix = substr(name, cut + 6)
    colon = 0
    for (i = length(prefix); i > 0; i--) if (substr(prefix, i, 1) == ":") { colon = i; break }
    mi = substr(prefix, colon + 1)
    if (colon == 0 || substr(prefix, 1, colon - 1) != id || mi !~ /^[0-9]+$/) {
      print "read name is not <id>:<MI>_CB:Z:<barcode>: " name > "/dev/stderr"; exit 1
    }
    tag = ""
    for (i = 12; i <= NF; i++) if (substr($i, 1, 5) == "CB:Z:") { tag = substr($i, 6); break }
    if (tag != suffix) { print "CB tag \"" tag "\" does not match read name " name > "/dev/stderr"; exit 1 }
    if (!(tag in seen)) { seen[tag] = 1; print tag > out }
  }
' || { log_error "read name / CB tag check failed for $par_id"; exit 1; }
sort -u -o "$out_barcodes" "$out_barcodes"
cell_barcodes="$meta_temp_dir/cell_barcodes.txt"
cut -f1 "$par_cell_list" | sort -u > "$cell_barcodes"
not_cells=$(comm -23 "$out_barcodes" "$cell_barcodes")
if [[ -n "$not_cells" ]]; then
  log_error "barcodes not in the cell list: $(head -5 <<< "$not_cells" | tr '\n' ' ')"
  exit 1
fi
log "✓ Every output barcode is in the cell list"
n_barcodes=$(wc -l < "$out_barcodes")
[[ "$n_barcodes" -ge 2 ]] || { log_error "expected at least 2 distinct barcodes, got $n_barcodes"; exit 1; }
log "✓ Every read name is $par_id:<MI>_CB:Z:<barcode> with a matching CB tag ($n_barcodes barcodes)"

in_barcodes="$meta_temp_dir/in_barcodes.txt"
samtools view "$par_input_bam" | awk -F'\t' '
  {
    cut = 0; rest = $1
    while ((k = index(rest, "_CB:Z:")) > 0) { cut += k; rest = substr(rest, k + 1) }
    if (cut > 0) print substr($1, cut + 6)
  }
' | sort -u > "$in_barcodes"
missing=$(comm -23 "$out_barcodes" "$in_barcodes")
if [[ -n "$missing" ]]; then
  log_error "barcodes not found in the input BAM: $(head -5 <<< "$missing" | tr '\n' ' ')"
  exit 1
fi
log "✓ Every output barcode occurs in the input BAM"

print_test_summary "consensus_dedup output contract ($par_id)"
