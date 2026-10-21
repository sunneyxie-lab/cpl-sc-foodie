#!/bin/bash

## VIASH START
## VIASH END

source "$meta_resources_dir/test_helpers.sh"
setup_test_env

# Reads spec lines "barcode mi[,mi...]" from stdin; "-" stands for a record without MI.
make_bam_py="$meta_temp_dir/make_bam.py"
cat > "$make_bam_py" <<'EOF'
import sys
import pysam

header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": "chr1", "LN": 100000}]}
with pysam.AlignmentFile(sys.argv[1], "wb", header=header) as bam:
    n = 0
    for line in sys.stdin:
        barcode, mis = line.split()
        for mi in mis.split(","):
            rec = pysam.AlignedSegment(bam.header)
            rec.query_name = f"r{n}"
            rec.reference_id = 0
            rec.reference_start = n % 1000
            rec.cigarstring = "10M"
            rec.query_sequence = "ACGTACGTAC"
            rec.query_qualities = pysam.qualitystring_to_array("IIIIIIIIII")
            rec.set_tag("CB", barcode, "Z")
            if mi != "-":
                rec.set_tag("MI", mi, "Z")
            bam.write(rec)
            n += 1
EOF
make_bam() {
  python3 "$make_bam_py" "$1"
}

# Spec line giving a barcode n distinct MIs.
groups_spec() {
  local barcode=$1 n=$2
  printf '%s %s\n' "$barcode" "$(seq -s, 0 $((n - 1)))"
}

run_case() {
  local name=$1 bam=$2
  shift 2
  out="$meta_temp_dir/$name"
  mkdir -p "$out"
  "$meta_executable" --input "$bam" --id "sample_$name" \
    --output_cell_list "$out/cell_list.tsv" \
    --output_cell_counts "$out/cell_counts.tsv" \
    --output_cell_qc "$out/cell_qc.tsv" \
    --output_knee_plot "$out/knee_plot.png" \
    "$@"
  for f in cell_list.tsv cell_counts.tsv cell_qc.tsv knee_plot.png; do
    check_file_not_empty "$out/$f" "$name $f"
  done
  png_magic=$(head -c4 "$out/knee_plot.png" | od -An -tx1 | tr -d ' ')
  [[ "$png_magic" == "89504e47" ]] || { log_error "$name: knee plot is not a PNG: $png_magic"; exit 1; }
  check_file_line_count "$out/cell_qc.tsv" 2 "$name QC summary"
}

# Column of the QC data row, by header name.
qc_field() {
  awk -F'\t' -v col="$2" 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == col) c = i } NR == 2 { print $c }' "$1"
}

check_state() {
  local got
  got=$(qc_field "$out/cell_qc.tsv" state)
  [[ "$got" == "$1" ]] || { log_error "expected state $1, got $got"; exit 1; }
  log "✓ State is $1"
}

check_cells() {
  local expected=$1 got
  got=$(cut -f1 "$out/cell_list.tsv" | sort | tr '\n' ' ')
  [[ "$got" == "$expected" ]] || { log_error "expected cells '$expected', got '$got'"; exit 1; }
  log "✓ Cell list is $expected"
}

log "Case 1: counting distinct MI per barcode"
bam1="$meta_temp_dir/count.bam"
printf '%s\n' "AAAA 0,0,0,1,1,-" "CCCC -,-" | make_bam "$bam1"
run_case count "$bam1"
printf 'CCCC\t0\nAAAA\t2\n' > "$meta_temp_dir/expected_counts.tsv"
cmp "$out/cell_counts.tsv" "$meta_temp_dir/expected_counts.tsv" || { log_error "unexpected counts: $(cat "$out/cell_counts.tsv")"; exit 1; }
log "✓ Counts are distinct MIs, and a barcode without MI counts 0"

spec="$meta_temp_dir/knee.spec"
: > "$spec"
for i in $(seq 1 30); do groups_spec "LOW$(printf '%02d' "$i")" $((i % 2 + 1)) >> "$spec"; done
for i in $(seq 0 9); do groups_spec "HIGH$i" $((40 + 2 * i)) >> "$spec"; done
knee_bam="$meta_temp_dir/knee.bam"
make_bam "$knee_bam" < "$spec"
# On a step, kneedle puts the elbow on the last point before the jump, the top low barcode (index 24 after the drop).
knee_cells="HIGH0 HIGH1 HIGH2 HIGH3 HIGH4 HIGH5 HIGH6 HIGH7 HIGH8 HIGH9 LOW29 "
top3="HIGH7 HIGH8 HIGH9 "
small_args=(--drop_lowest 5 --max_cells 20 --fallback_cells 3)

log "Case 2: a clear knee"
run_case normal "$knee_bam" "${small_args[@]}"
check_state Normal
check_cells "$knee_cells"
[[ "$(qc_field "$out/cell_qc.tsv" knee_index)" == 24 ]] || { log_error "expected knee_index 24"; exit 1; }
log "✓ Knee index is 24"

log "Case 3: a knee with max_cells or more cells falls back"
run_case max_cells "$knee_bam" --drop_lowest 5 --max_cells 5 --fallback_cells 3
check_state Abnormal
check_cells "$top3"
[[ "$(qc_field "$out/cell_qc.tsv" knee_index)" =~ ^[0-9]+$ ]] || { log_error "knee_index is not numeric"; exit 1; }
log "✓ Knee index is still reported"

log "Case 4: drop_lowest leaves no barcodes"
run_case too_few "$knee_bam" --drop_lowest 1000 --max_cells 20 --fallback_cells 3
check_state Abnormal
check_cells "$top3"
[[ "$(qc_field "$out/cell_qc.tsv" knee_index) $(qc_field "$out/cell_qc.tsv" knee_count)" == "NA NA" ]] ||
  { log_error "expected NA knee index and count"; exit 1; }
log "✓ Knee index and count are NA"

log "Case 5: fewer barcodes than fallback_cells"
bam5="$meta_temp_dir/two.bam"
{ groups_spec AAAA 3; groups_spec CCCC 5; } | make_bam "$bam5"
run_case two "$bam5"
check_state Abnormal
check_cells "AAAA CCCC "

print_test_summary "call_cells tests"
