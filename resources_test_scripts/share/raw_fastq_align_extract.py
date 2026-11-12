#!/usr/bin/env python3
"""Select BAM identities and preserve their original paired raw FASTQ records."""

import contextlib
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

from fixture_integrity import (
    CHROM, START, END, PREDICATE, SAMPLE_IDS,
    REFERENCE_LENGTH, REFERENCE_MD5, FixtureError,
    sha256_file, write_json, resources_root, staging_directory,
    verify_reference, common_cell_lists, verify_fixture_set, write_completion,
)

BAM_ROOT = Path("/lustre/grp/xxllab/fangj/03.project/02.snFOODIE/H0O260120CPL0071F1S/reRun/Align_bam")
FASTQ_ROOT = Path("/lustre/grp/xxllab/fangj/05.rawdata/02.invivo/CPL/147.H0O260120CPL0071F1S/Merge")
WHITELIST_SHA256 = "ecc4d583077e5d3b12d0bef70059820bf2579598c4177e8b1f926ca7e9ffee37"
COMMENT = re.compile(r"([12]):[YN]:[0-9]+:\S+\Z")
SELECTED_NAME = re.compile(r"(.+)_([12]):[YN]:[0-9]+:\S+\Z")


def normalize_selected_name(qname):
    prefix, marker, barcode = qname.rpartition("_CB:Z:")
    match = SELECTED_NAME.fullmatch(prefix)
    if not marker or not barcode or re.search(r"\s", barcode) or not match:
        raise FixtureError(f"unexpected selected QNAME grammar: {qname!r}")
    identity = match.group(1)
    if identity.startswith("@") or re.search(r"\s", identity):
        raise FixtureError(f"unexpected selected identity: {identity!r}")
    return identity, barcode


def within_interval(reference, position, mate_position):
    return reference == CHROM and START <= position <= END and START <= mate_position <= END


def select_sam_records(lines):
    """Use SAM identity, flags and coordinates only, never SAM sequence/quality."""
    selected = {}
    by_name = {}
    count = 0
    for number, line in enumerate(lines, 1):
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) < 11:
            raise FixtureError(f"selected SAM record {number}: expected at least 11 fields")
        qname = fields[0]
        try:
            flag, position, mate_position = int(fields[1]), int(fields[3]), int(fields[7])
        except ValueError as exc:
            raise FixtureError(f"selected SAM record {number}: invalid numeric field") from exc
        if not within_interval(fields[2], position, mate_position):
            raise FixtureError(f"selected SAM record {number}: outside inclusive predicate")
        mate = 1 if flag & 64 else 2
        if not flag & 1 or flag & (4 | 8) or bool(flag & 64) == bool(flag & 128):
            raise FixtureError(f"{qname}: expected paired, mapped first/second mate flags")
        identity, barcode = normalize_selected_name(qname)
        if identity in selected and selected[identity]["qname"] != qname:
            raise FixtureError(f"selected QNAME collision for identity {identity!r}")
        selected.setdefault(identity, {"qname": qname, "identity": identity, "barcode": barcode})
        records = by_name.setdefault(qname, {})
        if mate in records:
            raise FixtureError(f"{qname}: duplicate/extra selected mate {mate}")
        mate_reference = fields[2] if fields[6] == "=" else fields[6]
        records[mate] = (fields[2], position, mate_reference, mate_position)
        count += 1
    if not selected:
        raise FixtureError("BAM selection contains no pairs")
    for qname, records in by_name.items():
        if set(records) != {1, 2}:
            raise FixtureError(f"{qname}: missing selected BAM mate")
        first, second = records[1], records[2]
        if first[:2] != second[2:] or first[2:] != second[:2]:
            raise FixtureError(f"{qname}: selected BAM mate positions/references are not reciprocal")
    if count != 2 * len(selected):
        raise FixtureError("selected BAM record/pair counts disagree")
    return selected


def select_bam(bam):
    workers = max(0, int(os.environ.get("SLURM_CPUS_PER_TASK", "1")) - 1)
    command = ["samtools", "view", "--no-PG", "-@", str(workers), "-e", PREDICATE, str(bam)]
    with tempfile.TemporaryFile(mode="w+t") as errors:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=errors, text=True)
        try:
            selected = select_sam_records(process.stdout)
        except BaseException:
            process.stdout.close()
            process.terminate()
            process.wait()
            errors.seek(0)
            detail = errors.read().strip()
            if detail:
                print(detail, file=sys.stderr)
            raise
        finally:
            process.stdout.close()
        if process.wait():
            errors.seek(0)
            raise FixtureError(f"samtools selection failed for {bam}: {errors.read().strip()}")
    return selected


def read_record(handle, mate, number):
    header = handle.readline()
    if not header:
        return None
    lines = (header, handle.readline(), handle.readline(), handle.readline())
    if any(not line or not line.endswith(b"\n") for line in lines):
        raise FixtureError(f"R{mate} record {number}: incomplete four-line FASTQ")
    if not header.startswith(b"@") or not lines[2].startswith(b"+"):
        raise FixtureError(f"R{mate} record {number}: invalid header/separator")
    parts = header[1:].rstrip(b"\r\n").split()
    if len(parts) != 2:
        raise FixtureError(f"R{mate} record {number}: expected identity and Illumina comment")
    try:
        identity, comment = (part.decode("ascii") for part in parts)
    except UnicodeDecodeError as exc:
        raise FixtureError(f"R{mate} record {number}: non-ASCII header") from exc
    match = COMMENT.fullmatch(comment)
    if not identity or not match or int(match.group(1)) != mate:
        raise FixtureError(f"R{mate} record {number}: invalid mate-specific Illumina comment")
    sequence, quality = (line.rstrip(b"\r\n") for line in (lines[1], lines[3]))
    if not sequence or len(sequence) != len(quality):
        raise FixtureError(f"R{mate} record {number}: empty sequence or sequence/quality lengths differ")
    return identity, b"".join(lines)


def semantic_fields(record):
    lines = record.splitlines()
    return lines[0], lines[1], lines[3]


def paired_records(read1, read2, compressed=False):
    opener = gzip.open if compressed else open
    with opener(read1, "rb") as first, opener(read2, "rb") as second:
        number = 0
        while True:
            number += 1
            pair = [read_record(first, 1, number), read_record(second, 2, number)]
            if pair == [None, None]:
                return
            if None in pair:
                raise FixtureError(f"missing mate at record {number}")
            if pair[0][0] != pair[1][0]:
                raise FixtureError(f"desynchronized identities at record {number}")
            yield pair[0][0], (pair[0][1], pair[1][1])


def audit_raw_pairs(read1, read2, selected, audit1, audit2):
    """Scan original gzip sources through EOF and retain selected records verbatim."""
    seen = set()
    scanned = 0
    try:
        with open(audit1, "xb") as first, open(audit2, "xb") as second:
            for identity, pair in paired_records(read1, read2, compressed=True):
                scanned += 1
                if identity not in selected:
                    continue
                if identity in seen:
                    raise FixtureError(f"duplicate selected identity: {identity}")
                seen.add(identity)
                first.write(pair[0])
                second.write(pair[1])
    except (OSError, EOFError) as exc:
        raise FixtureError(f"original FASTQ audit failed for {read1}, {read2}: {exc}") from exc
    missing = selected.keys() - seen
    if missing:
        raise FixtureError(f"missing selected identities in original FASTQs: {len(missing)}")
    return {f"{field}_r{mate}": value
            for mate in (1, 2)
            for field, value in (("scanned", scanned), ("matched", len(seen)),
                                 ("missing", 0), ("duplicate", 0), ("written", len(seen)))}


def materialize_original_records(records, output1, output2):
    """Write validated original audit pairs in their original insertion order."""
    with gzip.open(output1, "xb") as first, gzip.open(output2, "xb") as second:
        for pair in records.values():
            for output, record in zip((first, second), pair):
                output.write(record)


def extract_raw_pairs(read1, read2, selected, output1, output2, *, cache_dir=None):
    """Select through CBQ, proving semantic identity against original-source audits."""
    threads = max(1, int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
    original_sources = {f"read{mate}": source_metadata(path) for mate, path in enumerate((read1, read2), 1)}
    if cache_dir is None:
        cache_dir = resources_root(Path(__file__).resolve().parents[2]) / ".raw_fastq_align_cbq"
    cache = Path(os.environ.get("RAW_FASTQ_ALIGN_CBQ_CACHE", cache_dir))
    cache.mkdir(parents=True, exist_ok=True)
    with contextlib.nullcontext(tempfile.mkdtemp(prefix=Path(read1).name + ".", dir=cache)) as temporary:
        work = Path(temporary)
        print(f"Retaining private CBQ workspace: {work}", flush=True)
        patterns = work / "patterns.tsv"
        patterns.write_text("".join(f"{index}\t{identity}\n" for index, identity in enumerate(selected)), encoding="ascii")
        audits = [work / f"original_R{mate}.fq" for mate in (1, 2)]
        counts_file = work / "counts.json"
        cbq = work / "source.cbq"
        print(f"Auditing original paired sources: {read1}, {read2}", flush=True)
        counts = audit_raw_pairs(read1, read2, selected, *audits)
        write_json(counts_file, counts)
        encode = ["bqtools", "encode", str(read1), str(read2), "--mode", "cbq",
                  "-T", str(threads), "-o", str(cbq)]
        print(f"Encoding original paired gzip sources to CBQ: {read1}, {read2}", flush=True)
        subprocess.run(encode, check=True)
        info_command = ["bqtools", "info", str(cbq), "--json"]
        info = json.loads(subprocess.check_output(info_command, text=True))
        if len(info) != 1 or any(info[0].get(field) != value for field, value in
                                  (("format", "CBQ"), ("paired", True), ("quality", True), ("headers", True),
                                   ("num_records", counts["scanned_r1"]))):
            raise FixtureError(f"CBQ metadata disagrees with complete source validation: {info}")
        if original_sources != {f"read{mate}": source_metadata(path) for mate, path in enumerate((read1, read2), 1)}:
            raise FixtureError("original FASTQ source changed during CBQ encoding")
        write_json(work / "cbq.completed.json", {
            "schema_version": 1, "completed": True, "format": "CBQ", "bqtools_version": "0.6.1",
            "sources": original_sources,
            "cbq": source_metadata(cbq), "metadata": info[0], "validated_source_counts": counts,
        })
        original = {}
        original_records = {}
        source_record_digests = [hashlib.sha256(), hashlib.sha256()]
        for identity, pair in paired_records(*audits):
            if identity not in selected or identity in original:
                raise FixtureError(f"invalid or duplicate original audit identity: {identity}")
            original[identity] = tuple(semantic_fields(record) for record in pair)
            original_records[identity] = pair
            for digest, record in zip(source_record_digests, pair):
                digest.update(record)
        if original.keys() != selected.keys():
            raise FixtureError("original audits do not cover selected identities")
        prefix = work / "candidates"
        grep = ["bqtools", "grep", str(cbq), "--header", "--fixed", "--sfile", str(patterns),
                "--format", "q", "--prefix", str(prefix), "-T", str(threads)]
        print(f"Extracting {len(selected)} selected original identities with bqtools grep", flush=True)
        subprocess.run(grep, check=True)
        seen = set()
        for identity, pair in paired_records(work / "candidates_R1.fq", work / "candidates_R2.fq"):
            # Fixed-string grep is a fast candidate search; equality defines selection.
            if identity not in selected:
                continue
            if identity in seen:
                raise FixtureError(f"duplicate selected identity in CBQ output: {identity}")
            seen.add(identity)
            fields = tuple(semantic_fields(record) for record in pair)
            if fields != original[identity]:
                raise FixtureError(f"CBQ changed an original header/comment, sequence or quality: {identity}")
        if seen != selected.keys():
            raise FixtureError(f"missing selected identities in CBQ output: {len(selected.keys() - seen)}")
        materialize_original_records(original_records, output1, output2)
        return {f"read{mate}": digest.hexdigest() for mate, digest in enumerate(source_record_digests, 1)}


def verify_output_records(output1, output2, expected_pairs, expected_digests):
    digests = [hashlib.sha256(), hashlib.sha256()]
    seen = set()
    with gzip.open(output1, "rb") as first, gzip.open(output2, "rb") as second:
        number = 0
        while True:
            pair = [read_record(first, 1, number + 1), read_record(second, 2, number + 1)]
            if pair == [None, None]:
                break
            if None in pair or pair[0][0] != pair[1][0] or pair[0][0] in seen:
                raise FixtureError("decoded output identities are incomplete, duplicated or desynchronized")
            seen.add(pair[0][0])
            for record, digest in zip(pair, digests):
                digest.update(record[1])
            number += 1
    actual = {f"read{mate}": digest.hexdigest() for mate, digest in enumerate(digests, 1)}
    if number != expected_pairs or actual != expected_digests:
        raise FixtureError("decoded output counts/digests differ from original selected records")


def source_metadata(path):
    path = Path(path).resolve()
    stat = path.stat()
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def tool_versions():
    output = subprocess.check_output(["samtools", "--version"], text=True)
    version = output.splitlines()[0]
    match = re.match(r"samtools ([0-9]+)\.([0-9]+)", version)
    if not match or tuple(map(int, match.groups())) < (1, 12):
        raise FixtureError(f"samtools >=1.12 required, found {version!r}")
    bqtools = subprocess.check_output(["bqtools", "--version"], text=True).strip()
    if bqtools != "bqtools 0.6.1":
        raise FixtureError(f"bqtools 0.6.1 required, found {bqtools!r}")


def verify_whitelist(path):
    digest = sha256_file(path)
    with open(path, encoding="ascii") as handle:
        entries = handle.read().splitlines()
    if digest != WHITELIST_SHA256 or len(entries) != 96 or any(len(entry) != 11 for entry in entries):
        raise FixtureError(f"whitelist identity differs from the fixed 96-entry source: {path}")
    return digest


def generate_sample(sample_id, bam, read1, read2, whitelist, output_dir, *, reference, cache_dir=None):
    if not sample_id or sample_id.startswith("@") or re.search(r"\s", sample_id):
        raise FixtureError(f"invalid sample ID: {sample_id!r}")
    # These fixed sample IDs also keep output names inside the staging directory.
    if sample_id not in SAMPLE_IDS:
        raise FixtureError(f"expected one of the two fixed sample IDs, got {sample_id!r}")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    stage = "preflight"
    try:
        if not isinstance(reference, dict) or tuple(reference.get(field) for field in ("name", "length", "md5")) != (CHROM, REFERENCE_LENGTH, REFERENCE_MD5):
            raise FixtureError("a validated full-chr14 reference context is required")
        verify_whitelist(whitelist)
        sources = {name: source_metadata(path) for name, path in (("bam", bam), ("read1", read1), ("read2", read2))}
        stage = "BAM identity selection"
        selected = select_bam(bam)
        stage = "original raw extraction"
        out1, out2 = [output_dir / f"{sample_id}.R{mate}.fastq.gz" for mate in (1, 2)]
        source_digests = extract_raw_pairs(read1, read2, selected, out1, out2, cache_dir=cache_dir)
        stage = "decoded output verification"
        verify_output_records(out1, out2, len(selected), source_digests)
        for name, path in (("bam", bam), ("read1", read1), ("read2", read2)):
            if source_metadata(path) != sources[name]:
                raise FixtureError(f"source changed during generation: {path}")
        mapping = output_dir / f"{sample_id}.identities.tsv"
        with open(mapping, "x", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["qname", "identity", "barcode"])
            writer.writerows((entry["qname"], entry["identity"], entry["barcode"]) for entry in selected.values())
    except Exception as exc:
        raise FixtureError(f"sample={sample_id} source_bam={bam} source_R1={read1} source_R2={read2} stage={stage}: {exc}") from exc


def generate_fixture_set(repo, resources_dir):
    repo = Path(repo).resolve()
    resources_dir = resources_root(repo, resources_dir)
    destination = resources_dir / "raw_fastq_align"
    if os.path.lexists(destination):
        raise FixtureError(f"refusing existing final destination: {destination}")
    common_manifest = verify_fixture_set("common", resources_dir)
    reference = verify_reference(resources_dir / "common/genome")
    common_cell_lists(resources_dir)
    tool_versions()
    whitelist = repo / "assets/barcodes.txt"
    verify_whitelist(whitelist)
    sources = [(sample_id, BAM_ROOT / f"{sample_id}.bam", FASTQ_ROOT / f"{sample_id}_ALL.R1.fastq.gz", FASTQ_ROOT / f"{sample_id}_ALL.R2.fastq.gz") for sample_id in SAMPLE_IDS]
    for sample_id, *paths in sources:
        for path in paths:
            if not path.is_file() or not os.access(path, os.R_OK):
                raise FixtureError(f"sample={sample_id} stage=preflight: cannot read source {path}")
    cache_dir = resources_dir / ".raw_fastq_align_cbq"
    with staging_directory(destination) as stage:
        for sample_id, bam, read1, read2 in sources:
            generate_sample(sample_id, bam, read1, read2, whitelist, stage, reference=reference, cache_dir=cache_dir)
        with open(stage / "samples.csv", "x", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(["id", "read1", "read2", "cell_list"])
            for sample_id in SAMPLE_IDS:
                writer.writerow([sample_id, f"{sample_id}.R1.fastq.gz", f"{sample_id}.R2.fastq.gz", f"../common/cell_lists/{sample_id}.cell_list.tsv"])
        write_completion("raw_fastq_align", stage, common_manifest=common_manifest)
        verify_fixture_set("raw_fastq_align", resources_dir, directory=stage)
    return destination
