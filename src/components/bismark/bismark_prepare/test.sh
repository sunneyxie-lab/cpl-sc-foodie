#!/bin/bash

## VIASH START
## VIASH END

source "$meta_resources_dir/test_helpers.sh"
setup_test_env

fasta="$meta_temp_dir/chr1.fa"
awk 'BEGIN {
  srand(1)
  print ">chr1"
  for (i = 0; i < 5000; i++) {
    printf "%s", substr("ACGT", int(rand() * 4) + 1, 1)
    if (i % 60 == 59) printf "\n"
  }
  printf "\n"
}' > "$fasta"

log "Preparing the Bismark genome"
out_dir="$meta_temp_dir/genome"
"$meta_executable" --input "$fasta" --output "$out_dir"

check_dir_exists "$out_dir" "genome directory"
check_file_exists "$out_dir/chr1.fa" "copied FASTA"
[[ ! -L "$out_dir/chr1.fa" ]] || { log_error "FASTA is a symlink, expected a copy"; exit 1; }
for conv in CT_conversion GA_conversion; do
  check_dir_exists "$out_dir/Bisulfite_Genome/$conv" "$conv index directory"
  n_bt2=$(find "$out_dir/Bisulfite_Genome/$conv" -name '*.bt2' | wc -l)
  [[ "$n_bt2" -gt 0 ]] || { log_error "no .bt2 files in $conv"; exit 1; }
  log "✓ $conv has $n_bt2 Bowtie 2 index files"
done

print_test_summary "bismark_prepare tests"
