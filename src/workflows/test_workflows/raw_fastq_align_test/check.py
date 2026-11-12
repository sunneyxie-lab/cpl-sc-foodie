#!/usr/bin/env python3
"""Content checks for the test-only raw FASTQ producer contract."""

import argparse
import csv
import gzip
import html
import math
from pathlib import Path
import re
import subprocess
import sys
import zipfile


class CheckError(Exception):
    pass


def require(condition, message):
    if not condition:
        raise CheckError(message)


def raw_records(path, mate):
    with gzip.open(path, "rb") as handle:
        number = 0
        while True:
            header = handle.readline()
            if not header:
                return
            number += 1
            lines = [header] + [handle.readline() for _ in range(3)]
            label = f"raw R{mate} record {number}"
            require(all(line.endswith(b"\n") for line in lines), f"{label}: incomplete four-line record")
            head, sequence, plus, quality = [line.rstrip(b"\r\n") for line in lines]
            require(head.startswith(b"@") and plus.startswith(b"+"), f"{label}: malformed FASTQ markers")
            fields = head[1:].split()
            require(len(fields) == 2 and re.fullmatch(rb"[12]:[YN]:[0-9]+:\S+", fields[1]),
                    f"{label}: missing or malformed Illumina mate comment")
            require(fields[1].startswith(str(mate).encode() + b":"), f"{label}: wrong mate comment")
            require(sequence and len(sequence) == len(quality), f"{label}: sequence/quality length mismatch")
            require(all(33 <= value <= 126 for value in quality), f"{label}: invalid Phred+33 quality")
            require(head.isascii() and sequence.isascii(), f"{label}: non-ASCII header or sequence")
            yield fields[0].decode("ascii"), sequence, quality


def inspect_raw(read1, read2):
    from itertools import zip_longest

    identities = set()
    metrics = [{"num_seqs": 0, "sum_len": 0, "min_len": None, "max_len": 0,
                "gc": 0, "q20": 0, "q30": 0} for _ in range(2)]
    for first, second in zip_longest(raw_records(read1, 1), raw_records(read2, 2)):
        require(first is not None and second is not None, "raw inputs: unequal mate counts")
        require(first[0] == second[0], f"raw inputs: desynchronized identity {first[0]} / {second[0]}")
        require(first[0] not in identities, f"raw inputs: duplicate identity {first[0]}")
        identities.add(first[0])
        for record, metric in zip((first, second), metrics):
            _, sequence, quality = record
            size = len(sequence)
            metric["num_seqs"] += 1
            metric["sum_len"] += size
            metric["min_len"] = size if metric["min_len"] is None else min(metric["min_len"], size)
            metric["max_len"] = max(metric["max_len"], size)
            metric["gc"] += sum(base in b"GCgc" for base in sequence)
            metric["q20"] += sum(value >= 53 for value in quality)
            metric["q30"] += sum(value >= 63 for value in quality)
    require(identities, "raw inputs: no paired records")
    return identities, metrics


def number(value, label):
    try:
        result = float(value.replace(",", ""))
    except (ValueError, AttributeError) as error:
        raise CheckError(f"{label}: not a numeric metric: {value!r}") from error
    require(math.isfinite(result) and result >= 0, f"{label}: invalid nonnegative metric {value!r}")
    return result


def stats_table(path, expected_rows, label):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        required = {"file", "format", "type", "num_seqs", "sum_len", "min_len", "avg_len", "max_len",
                    "Q20(%)", "Q30(%)", "GC(%)", "N50"}
        require(required <= set(fields) and len(fields) == len(set(fields)), f"{label}: missing/duplicate full-statistics header")
        rows = list(reader)
    require(len(rows) == expected_rows, f"{label}: expected {expected_rows} data rows, got {len(rows)}")
    for index, row in enumerate(rows, 1):
        require(None not in row and all(value is not None and value.strip() for value in row.values()),
                f"{label} row {index}: malformed or empty table field")
        require(row["format"] == "FASTQ" and row["type"] == "DNA", f"{label} row {index}: not DNA FASTQ statistics")
        for field in required - {"file", "format", "type"}:
            number(row[field], f"{label} row {index} {field}")
        for field in ("num_seqs", "sum_len", "min_len", "max_len", "N50"):
            require(number(row[field], label).is_integer(), f"{label} row {index}: {field} is not integral")
        require(number(row["num_seqs"], label) > 0, f"{label} row {index}: empty FASTQ")
        require(number(row["min_len"], label) <= number(row["avg_len"], label) <= number(row["max_len"], label),
                f"{label} row {index}: inconsistent length metrics")
        for field in ("Q20(%)", "Q30(%)", "GC(%)"):
            require(number(row[field], label) <= 100, f"{label} row {index}: {field} exceeds 100")
        require(number(row["Q30(%)"], label) <= number(row["Q20(%)"], label),
                f"{label} row {index}: Q30 exceeds Q20")
    return rows


def match_raw_stats(rows, reads, metrics):
    for mate, (row, read, metric) in enumerate(zip(rows, reads, metrics), 1):
        label = f"raw statistics R{mate}"
        require(Path(row["file"]).name == Path(read).name, f"{label}: input association/order mismatch")
        for key in ("num_seqs", "sum_len", "min_len", "max_len"):
            require(number(row[key], label) == metric[key], f"{label}: {key} differs from decoded raw content")
        expected = {"avg_len": metric["sum_len"] / metric["num_seqs"],
                    "GC(%)": 100 * metric["gc"] / metric["sum_len"],
                    "Q20(%)": 100 * metric["q20"] / metric["sum_len"],
                    "Q30(%)": 100 * metric["q30"] / metric["sum_len"]}
        for key, value in expected.items():
            # SeqKit rounds quality percentages to two decimals before integer TSV formatting.
            if key in ("Q20(%)", "Q30(%)"):
                value = round(value, 2)
            decimals = len(row[key].partition(".")[2])
            tolerance = 0.5 * 10 ** -decimals + 1e-9
            require(abs(number(row[key], label) - value) <= tolerance, f"{label}: {key} differs from decoded raw content")


def inspect_qc(html_path, zip_path, row, label, native_filename=None):
    text = html.unescape(Path(html_path).read_text())
    require(re.search(r"<html\b", text, re.I) and re.search(r"<body\b", text, re.I) and "FastQC" in text,
            f"{label}: malformed FastQC HTML")
    accepted_files = {Path(row["file"]).name}
    if native_filename:
        accepted_files.add(native_filename)
    with zipfile.ZipFile(zip_path) as archive:
        require(archive.testzip() is None, f"{label}: corrupt ZIP member")
        data_files = [name for name in archive.namelist() if name.endswith("/fastqc_data.txt")]
        summaries = [name for name in archive.namelist() if name.endswith("/summary.txt")]
        require(len(data_files) == len(summaries) == 1, f"{label}: expected one FastQC data/summary pair")
        require(Path(data_files[0]).parent == Path(summaries[0]).parent, f"{label}: mixed FastQC ZIP report roots")
        data = archive.read(data_files[0]).decode("utf-8")
        summary = archive.read(summaries[0]).decode("utf-8")
    require(data.startswith("##FastQC\t"), f"{label}: malformed FastQC data header")
    match = re.search(r"^>>Basic Statistics\t(?:pass|warn|fail)\n(.*?)^>>END_MODULE", data, re.M | re.S)
    require(match is not None, f"{label}: missing Basic Statistics module")
    basic = dict(line.split("\t", 1) for line in match[1].splitlines() if "\t" in line and not line.startswith("#"))
    report_file = basic.get("Filename")
    require(report_file in accepted_files, f"{label}: ZIP input filename association/order mismatch")
    require(report_file in text, f"{label}: HTML input filename association mismatch")
    require(number(basic.get("Total Sequences", ""), label) == number(row["num_seqs"], label),
            f"{label}: FastQC sequence count differs from statistics")
    lengths = basic.get("Sequence length", "").split("-")
    require(len(lengths) in (1, 2) and number(lengths[0], label) == number(row["min_len"], label)
            and number(lengths[-1], label) == number(row["max_len"], label),
            f"{label}: FastQC sequence lengths differ from statistics")
    require(abs(number(basic.get("%GC", ""), label) - number(row["GC(%)"], label)) <= 1.01,
            f"{label}: FastQC GC differs from statistics")
    summary_rows = [line.split("\t") for line in summary.splitlines() if line]
    require(summary_rows and all(len(parts) == 3 and parts[0] in ("PASS", "WARN", "FAIL")
                                 and parts[2] == report_file for parts in summary_rows),
            f"{label}: malformed or misassociated FastQC summary")


def inspect_bam(path, raw_identities, trimmed_pairs):
    quick = subprocess.run(["samtools", "quickcheck", str(path)], capture_output=True, text=True)
    require(quick.returncode == 0, "BAM: samtools quickcheck failed")
    pairs = {}
    process = subprocess.Popen(["samtools", "view", str(path)], stdout=subprocess.PIPE, text=True)
    try:
        for line_number, line in enumerate(process.stdout, 1):
            fields = line.rstrip("\n").split("\t")
            require(len(fields) >= 11, f"BAM record {line_number}: malformed SAM record")
            name = fields[0]
            prefix, marker, barcode = name.rpartition("_CB:Z:")
            require(marker and prefix and barcode and not re.search(r"\s", barcode),
                    f"BAM record {line_number}: missing terminal nonempty barcode suffix: {name}")
            comment = re.fullmatch(r"(.+)_[12]:[YN]:[0-9]+:\S+", prefix)
            identity = comment.group(1) if comment else prefix
            require(identity in raw_identities, f"BAM record {line_number}: identity not present in this sample's raw inputs: {name}")
            require(not any(field.startswith("CB:") for field in fields[11:]),
                    f"BAM record {line_number}: producer added a CB auxiliary field: {name}")
            flag = int(fields[1])
            require(flag & 1 and not flag & 12 and bool(flag & 64) != bool(flag & 128),
                    f"BAM record {line_number}: not a mapped paired first/second mate: {name}")
            mate = 1 if flag & 64 else 2
            entry = pairs.setdefault(identity, {})
            require(mate not in entry, f"BAM: duplicate mate {mate} for identity {identity}")
            entry[mate] = (name, barcode, fields[2], int(fields[3]), fields[6], int(fields[7]))
        returncode = process.wait()
        require(returncode == 0, f"BAM: samtools view failed with exit {returncode}; see stderr above")
    finally:
        if process.poll() is None:
            process.kill()
            process.wait()
        process.stdout.close()
    require(pairs, "BAM: no aligned paired identities")
    require(len(pairs) <= trimmed_pairs, "BAM: more aligned identities than trimmed input pairs")
    for identity, mates in pairs.items():
        require(set(mates) == {1, 2}, f"BAM: missing mate for identity {identity}")
        first, second = mates[1], mates[2]
        require(first[:2] == second[:2], f"BAM: mate QNAME/barcode mismatch for identity {identity}")
        for record, other in ((first, second), (second, first)):
            reference = record[2] if record[4] == "=" else record[4]
            require(reference == other[2] and record[5] == other[3], f"BAM: nonreciprocal mate positions for identity {identity}")
    return len(pairs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    names = ("id", "read1", "read2", "output_bam", "output_bismark_report",
             "output_raw_stats", "output_extracted_trimmed_stats")
    for name in names:
        parser.add_argument("--" + name, required=True)
    for stage in ("raw", "extracted", "trimmed"):
        for extension in ("html", "zip"):
            parser.add_argument(f"--output_{stage}_fastqc_{extension}", required=True)
    args = parser.parse_args()
    require(re.fullmatch(r"[^@\s]\S*", args.id), "sample ID: invalid whitespace or leading @")
    reads = [args.read1, args.read2]
    identities, metrics = inspect_raw(*reads)
    raw = stats_table(args.output_raw_stats, 2, "raw statistics")
    processed = stats_table(args.output_extracted_trimmed_stats, 4, "extracted/trimmed statistics")
    match_raw_stats(raw, reads, metrics)
    for mate in (1, 2):
        require(Path(processed[mate - 1]["file"]).name == f"extracted_R{mate}.fastq.gz",
                f"extracted statistics R{mate}: filename association/order mismatch")
        trimmed_name = Path(processed[mate + 1]["file"]).name
        require(re.fullmatch(rf"trimmed(?:_R{mate}\.fastq|_val_{mate}\.fq)\.gz", trimmed_name),
                f"trimmed statistics R{mate}: filename association/order mismatch")
    for stage, rows in (("raw", raw), ("extracted", processed[:2]), ("trimmed", processed[2:])):
        require(number(rows[0]["num_seqs"], stage) == number(rows[1]["num_seqs"], stage),
                f"{stage} statistics: unequal mate counts")
        groups = []
        for extension in ("html", "zip"):
            paths = getattr(args, f"output_{stage}_fastqc_{extension}").split(";")
            require(len(paths) == 2 and all(paths) and paths[0] != paths[1],
                    f"{stage} FastQC {extension}: expected two distinct files in R1/R2 order")
            groups.append(paths)
        for mate, (html_path, zip_path, row) in enumerate(zip(*groups, rows), 1):
            native_filename = f"trimmed_val_{mate}.fq.gz" if stage == "trimmed" else None
            inspect_qc(html_path, zip_path, row, f"{stage} FastQC R{mate}", native_filename)
    raw_pairs, extracted_pairs, trimmed_pairs = [int(number(rows[0]["num_seqs"], "pair count"))
                                                for rows in (raw, processed[:2], processed[2:])]
    require(0 < trimmed_pairs <= extracted_pairs <= raw_pairs, "statistics: pair counts increase during extraction/trimming")
    for mate in range(2):
        require(number(processed[mate + 2]["sum_len"], "trimmed length") <= number(processed[mate]["sum_len"], "extracted length"),
                f"statistics: trimmed R{mate + 1} has more bases than extracted input")
    report = Path(args.output_bismark_report).read_text()
    totals = re.findall(r"^\s*Sequence pairs analy[sz]ed in total:\s*([0-9][0-9,]*)\s*$", report, re.M)
    require(len(totals) == 1, "Bismark report: missing or ambiguous analysed pair total")
    require(int(totals[0].replace(",", "")) == trimmed_pairs, "Bismark report: analysed pair total differs from trimmed statistics")
    aligned_pairs = inspect_bam(args.output_bam, identities, trimmed_pairs)
    unique = re.findall(r"^\s*Number of paired-end alignments with a unique best hit:\s*([0-9][0-9,]*)\s*$", report, re.M)
    require(len(unique) == 1, "Bismark report: missing or ambiguous unique paired alignment count")
    require(int(unique[0].replace(",", "")) == aligned_pairs, "Bismark report: unique paired alignment count differs from BAM")
    print(f"PASS raw_fastq_align ({args.id}): raw={raw_pairs}, extracted={extracted_pairs}, trimmed={trimmed_pairs}, aligned={aligned_pairs}; all ten output fields validated")


if __name__ == "__main__":
    try:
        main()
    except (CheckError, OSError, ValueError, KeyError, zipfile.BadZipFile, EOFError) as error:
        print(f"FAIL raw_fastq_align: {error}", file=sys.stderr)
        sys.exit(1)
