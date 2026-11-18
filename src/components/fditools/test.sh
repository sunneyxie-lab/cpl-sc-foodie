#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail
source "$meta_resources_dir/test_helpers.sh"
setup_test_env

python3 - "$meta_temp_dir" <<'PY'
import gzip
import pathlib
import sys

temp = pathlib.Path(sys.argv[1])
root = temp / "synthetic reads"
root.mkdir()
barcodes = ("ACGTACGTACG", "TGCATGCATGC")
(root / "whitelist.txt").write_text("\n".join(barcodes) + "\n")
# Interleave rejected pairs with exact and corrected matches to check pair order.
cases = [
    ("exact_ab", barcodes[0], barcodes[1], barcodes[0] + barcodes[1]),
    ("reject_r1", "C" * 11, barcodes[0], None),
    ("correct_r1", "TCGTACGTACG", barcodes[1], barcodes[0] + barcodes[1]),
    ("exact_ba", barcodes[1], barcodes[0], barcodes[1] + barcodes[0]),
    ("reject_r2", barcodes[1], "C" * 11, None),
    ("correct_r2", barcodes[0], "AGCATGCATGC", barcodes[0] + barcodes[1]),
]
for mate in (1, 2):
    expected = []
    with gzip.open(root / f"R{mate}.fastq.gz", "wt") as handle:
        for index, (name, bc1, bc2, barcode) in enumerate(cases):
            payload = ("GATTACA" if mate == 1 else "CTAGCTA") * (index + 1)
            quality = "".join(chr(40 + i % 30) for i in range(len(payload)))
            header = f"@{name} {mate}:N:0:INDEX"
            sequence = (bc1 if mate == 1 else bc2) + "T" * 19 + payload
            handle.write(f"{header}\n{sequence}\n+\n{'!' * 11}{'#' * 19}{quality}\n")
            if barcode is not None:
                expected.append(f"{header} CB:Z:{barcode}\n{payload}\n+\n{quality}\n")
    (root / f"expected.R{mate}.fastq").write_text("".join(expected))
PY

for cpus in 1 2 3 6; do
  prefix="$meta_temp_dir/synthetic.$cpus"
  input="$meta_temp_dir/synthetic reads"
  log "Checking synthetic extraction: $cpus CPUs"
  "$meta_executable" ---cpus "$cpus" \
    -a "$input/R1.fastq.gz" -b "$input/R2.fastq.gz" \
    -w "$input/whitelist.txt" \
    -c "$prefix.wrapper.R1.gz" -d "$prefix.wrapper.R2.gz"
  python3 - "$prefix" "$input" <<'PY'
import gzip
from pathlib import Path
import sys

prefix, root = sys.argv[1], Path(sys.argv[2])
for mate in (1, 2):
    with gzip.open(f"{prefix}.wrapper.R{mate}.gz", "rb") as handle:
        actual = handle.read()
    expected = (root / f"expected.R{mate}.fastq").read_bytes()
    assert actual == expected, f"R{mate}: unexpected extraction, correction, filtering, or record order"
PY
done
print_test_summary "fditools tests"
