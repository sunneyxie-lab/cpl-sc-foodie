#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail
source "$meta_resources_dir/test_helpers.sh"
setup_test_env

python3 - "$meta_resources_dir" "$meta_temp_dir" <<'PY'
import gzip
import importlib.util
import itertools
import json
import pathlib
import sys

spec = importlib.util.spec_from_file_location("selector", pathlib.Path(sys.argv[1]) / "select_legacy_adapter.py")
s = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s)
for counts, expected in [
    ((0, 0, 0), "illumina"), ((5, 1, 2), "illumina"),
    ((1, 5, 2), "nextera"), ((1, 2, 5), "small_rna"),
    ((5, 5, 1), "illumina"), ((5, 1, 5), "illumina"),
    ((1, 5, 5), "nextera"), ((5, 5, 5), "illumina"),
]:
    assert s.select_adapter(dict(zip(s.PROBES, counts))) == expected
probes = list(s.PROBES.values())
assert s.count_adapters([probes[0] * 3, probes[0].lower()]) == dict(illumina=1, nextera=0, small_rna=0)
assert s.select_adapter(s.count_adapters(["AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA"] * 50)) == "illumina"
consumed = 0

def records():
    global consumed
    for sequence in itertools.chain(itertools.repeat("A", 999_998), [probes[1], probes[2], probes[2]]):
        consumed += 1
        yield sequence

counts = s.count_adapters(records())
assert counts == dict(illumina=0, nextera=1, small_rna=0), counts
assert consumed == 1_000_000, consumed

root = pathlib.Path(sys.argv[2])
core = "ACCTATCACCTATCACCTATCACCTATCACCT"  # No adapter or poly-G motif.
cases = [
    ("adapter", core + probes[0], "I" * (len(core) + len(probes[0])), core),
    ("quality", core + "ATATATATAT", "I" * len(core) + "!" * 10, core),
    ("cutoff20", core[:20], "I" * 20, core[:20]),
    ("cutoff19", core[:19], "I" * 19, None),
]
# More than ten poly-G records trigger Rust's automatic poly-G detection.
cases += [(f"polyG{i}", core + "G" * 20, "I" * (len(core) + 20), core + "G" * 20) for i in range(12)]
expected = []
for name, sequence, quality, retained in cases:
    if retained is not None:
        expected.append([name, retained, "I" * len(retained)])
for mate in (1, 2):
    with gzip.open(root / f"R{mate}.fastq.gz", "wt") as handle:
        for name, sequence, quality, _ in cases:
            if name == "cutoff19" and mate == 2:
                sequence, quality = core, "I" * len(core)
            handle.write(f"@{name} {mate}:N:0:INDEX CB:Z:ACGTACGT\n{sequence}\n+\n{quality}\n")
(root / "expected.json").write_text(json.dumps(expected))
PY

log 'Rejecting malformed report templates without changing outputs'
python3 - "$meta_executable" "$meta_temp_dir" <<'PY'
import pathlib
import subprocess
import sys

executable, temp = sys.argv[1:]
root = pathlib.Path(temp)
for argument, suffix in (("fastqc_html", "html"), ("fastqc_zip", "zip")):
    for case, pattern in (
        ("no-wildcard", f"QC.{suffix}"),
        ("multiple-wildcards", f"QC.*.*.{suffix}"),
        ("multiple-values", f"QC.*.{suffix};other.*.{suffix}"),
    ):
        directory = root / f"invalid-{argument}-{case}"
        directory.mkdir()
        templates = {
            "fastqc_html": str(directory / "QC.*.html"),
            "fastqc_zip": str(directory / "QC.*.zip"),
        }
        templates[argument] = ";".join(str(directory / part) for part in pattern.split(";"))
        read1 = directory / "R1.fq.gz"
        read2 = directory / "R2.fq.gz"
        paths = {read1, read2}
        for template in templates.values():
            for part in template.split(";"):
                paths.update(pathlib.Path(part.replace("*", str(index))) for index in (0, 1))
        for path in paths:
            path.write_bytes(b"existing output\n")
        before = {path.relative_to(directory): path.read_bytes()
                  for path in directory.rglob("*") if path.is_file()}
        result = subprocess.run(
            [
                executable,
                "--input", f"{root / 'R1.fastq.gz'};{root / 'R2.fastq.gz'}",
                "--output_read1", str(read1),
                "--output_read2", str(read2),
                "--fastqc_html", templates["fastqc_html"],
                "--fastqc_zip", templates["fastqc_zip"],
            ],
            capture_output=True, text=True,
        )
        assert result.returncode != 0, (argument, case, result.stdout, result.stderr)
        after = {path.relative_to(directory): path.read_bytes()
                 for path in directory.rglob("*") if path.is_file()}
        assert after == before, (argument, case)
PY

log 'Checking real Rust trimming and bundled QC'
"$meta_executable" \
  --input "$meta_temp_dir/R1.fastq.gz;$meta_temp_dir/R2.fastq.gz" \
  --output_read1 "$meta_temp_dir/nested outputs/read 1/trimmed.fq.gz" \
  --output_read2 "$meta_temp_dir/nested outputs/read 2/trimmed.fq.gz" \
  --fastqc_html "$meta_temp_dir/nested outputs/html reports/mate */QC.html" \
  --fastqc_zip "$meta_temp_dir/nested outputs/zip reports/mate */QC.zip"

# Control: automatic poly-G must change the same reads when not disabled.
trim_galore --paired -q 20 --stringency 2 -e 0.1 --illumina \
  --basename control --output_dir "$meta_temp_dir/control" --cores 1 \
  "$meta_temp_dir/R1.fastq.gz" "$meta_temp_dir/R2.fastq.gz"

python3 - "$meta_temp_dir" <<'PY'
import gzip
import json
import pathlib
import sys
import zipfile

root = pathlib.Path(sys.argv[1])
expected = [tuple(record) for record in json.loads((root / "expected.json").read_text())]

def read(path):
    with path.open("rb") as raw:
        assert raw.read(2) == b"\x1f\x8b", path
    records = []
    with gzip.open(path, "rt") as handle:
        while header := handle.readline():
            sequence = handle.readline().rstrip()
            assert handle.readline().startswith("+")
            quality = handle.readline().rstrip()
            assert len(sequence) == len(quality)
            records.append((header.strip()[1:], sequence, quality))
    return records

for mate in (1, 2):
    mate_expected = [(f"{name} {mate}:N:0:INDEX CB:Z:ACGTACGT", seq, qual) for name, seq, qual in expected]
    observed = read(root / "nested outputs" / f"read {mate}" / "trimmed.fq.gz")
    assert observed == mate_expected, (mate, observed, mate_expected)
    control = dict((name, (seq, qual)) for name, seq, qual in read(root / "control" / f"control_val_{mate}.fq.gz"))
    for name, seq, qual in mate_expected:
        if name.startswith("polyG"):
            assert control[name] == (seq[:-20], qual[:-20]), control[name]
    html = (root / "nested outputs" / "html reports" / f"mate {mate - 1}" / "QC.html").read_text()
    assert "<html" in html.lower() and "fastqc" in html.lower()
    with zipfile.ZipFile(root / "nested outputs" / "zip reports" / f"mate {mate - 1}" / "QC.zip") as archive:
        assert archive.testzip() is None
        data = [name for name in archive.namelist() if name.endswith("/fastqc_data.txt")]
        assert len(data) == 1
        assert b"Total Sequences\t15" in archive.read(data[0])
        assert f"trimmed_val_{mate}.fq.gz".encode() in archive.read(data[0])
PY
print_test_summary 'trim_galore tests'
