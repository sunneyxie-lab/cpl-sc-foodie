#!/bin/bash

## VIASH START
## VIASH END

source "$meta_resources_dir/test_helpers.sh"
setup_test_env

genome_dir="$meta_temp_dir/genome"
mkdir -p "$genome_dir"
awk 'BEGIN {
  srand(1)
  print ">chr1"
  for (i = 0; i < 5000; i++) {
    printf "%s", substr("ACGT", int(rand() * 4) + 1, 1)
    if (i % 60 == 59) printf "\n"
  }
  printf "\n"
}' > "$genome_dir/chr1.fa"
bismark prepare "$genome_dir"

# R1 is the forward 50 bp at p, R2 the reverse complement of the 50 bp at p + 150.
n_pairs=10
r1="$meta_temp_dir/R1.fq"
r2="$meta_temp_dir/R2.fq"
awk -v n="$n_pairs" -v r1="$r1" -v r2="$r2" '
  function revcomp(s,   i, out) {
    out = ""
    for (i = length(s); i > 0; i--) out = out comp[substr(s, i, 1)]
    return out
  }
  !/^>/ { g = g $0 }
  END {
    comp["A"] = "T"; comp["C"] = "G"; comp["G"] = "C"; comp["T"] = "A"
    q = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
    for (i = 0; i < n; i++) {
      p = 201 + i * 400
      printf "@H0O_001:%d\n%s\n+\n%s\n", i, substr(g, p, 50), q > r1
      printf "@H0O_001:%d\n%s\n+\n%s\n", i, revcomp(substr(g, p + 150, 50)), q > r2
    }
  }
' "$genome_dir/chr1.fa"
gzip "$r1" "$r2"

log "Aligning $n_pairs read pairs"
out_bam="$meta_temp_dir/out.bam"
report="$meta_temp_dir/report.txt"
"$meta_executable" \
  --genome "$genome_dir" \
  --mates1 "$r1.gz" \
  --mates2 "$r2.gz" \
  --maxins 2000 \
  --output "$out_bam" \
  --report "$report"

check_file_exists "$out_bam" "Bismark BAM"
magic=$(head -c2 "$out_bam" | od -An -tx1)
[[ "$magic" == *"1f 8b"* ]] || { log_error "output is not BGZF: $magic"; exit 1; }
n_named=$(zcat "$out_bam" | grep -a -c 'H0O_001:' || true)
[[ "$n_named" -gt 0 ]] || { log_error "consensus read names were not kept"; exit 1; }
log "✓ BAM keeps the H0O_001:<i> read names"

check_file_exists "$report" "alignment report"
check_file_matches_regex "$report" "Sequence pairs analysed in total:[[:space:]]+$n_pairs\$" "alignment report"
check_file_contains "$report" "Mapping efficiency:" "alignment report"

print_test_summary "bismark_align tests"
