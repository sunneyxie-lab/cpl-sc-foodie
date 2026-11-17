#!/usr/bin/env python3
"""Run the independent Perl baseline and compare every retained paired record."""

import argparse
import gzip
import hashlib
from itertools import zip_longest
import json
from pathlib import Path
import platform
import re
import shlex
import shutil
import subprocess
import sys


PERL_SOURCE = "https://raw.githubusercontent.com/FelixKrueger/TrimGalore/0.6.11/trim_galore"
PROBES = {"illumina": "AGATCGGAAGAGC", "nextera": "CTGTCTCTTATA", "small_rna": "TGGAATTCTCGG"}


class CheckError(Exception):
    pass


def require(condition, message):
    if not condition:
        raise CheckError(message)


def checksum(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def record_digest(digest, record):
    # Length-safe canonical JSON tuples, independent of gzip/header comments/plus lines.
    digest.update((json.dumps(record, ensure_ascii=True, separators=(",", ":")) + "\n").encode("ascii"))


def records(path, label):
    with open(path, "rb") as handle:
        require(handle.read(2) == b"\x1f\x8b", f"{label}: not a gzip FASTQ")
    with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
        index = 0
        while True:
            header = handle.readline()
            if not header:
                return
            index += 1
            lines = [header] + [handle.readline() for _ in range(3)]
            require(all(line.endswith("\n") for line in lines), f"{label} record {index}: incomplete four-line FASTQ")
            head, sequence, plus, quality = [line.rstrip("\r\n") for line in lines]
            require(head.startswith("@") and plus.startswith("+"), f"{label} record {index}: malformed FASTQ markers")
            fields = head[1:].split()
            require(fields, f"{label} record {index}: empty identity")
            barcode = re.search(r"\sCB:Z:(\S+)$", head[1:])
            require(barcode is not None, f"{label} record {index}: missing terminal nonempty CB:Z: comment: {head}")
            name = (fields[0], barcode.group(1))
            require(sequence and len(sequence) == len(quality), f"{label} record {index}: sequence/quality length mismatch")
            require(all(33 <= ord(value) <= 126 for value in quality), f"{label} record {index}: invalid Phred+33 quality")
            yield name, sequence, quality, head[1:]


def pairs(read1, read2, label):
    seen = set()
    for index, (first, second) in enumerate(zip_longest(records(read1, f"{label} R1"), records(read2, f"{label} R2")), 1):
        require(first is not None and second is not None, f"{label}: unequal mate counts at pair {index}")
        require(first[0] == second[0], f"{label}: mate name/barcode mismatch at pair {index}: {first[0]} / {second[0]}")
        require(first[0] not in seen, f"{label}: duplicate identity {first[0]}")
        seen.add(first[0])
        yield first[0], first[1], first[2], second[1], second[2], first[3], second[3]


def inspect_inputs(read1, read2):
    identities = {}
    probe_counts = {name: 0 for name in PROBES}
    content = hashlib.sha256()
    mate_content = [hashlib.sha256(), hashlib.sha256()]
    for index, record in enumerate(pairs(read1, read2, "extracted inputs")):
        identities[record[0]] = {"position": index, "headers": record[5:7]}
        record_digest(content, record[:5])
        record_digest(mate_content[0], record[:3])
        record_digest(mate_content[1], (record[0], record[3], record[4]))
        if index < 999999:
            for name, probe in PROBES.items():
                probe_counts[name] += probe in record[1]
    require(identities, "extracted inputs: no pairs")
    # Explicit priority implements Illumina ties, then Nextera over small RNA.
    selected = max(PROBES, key=lambda name: (probe_counts[name], {"illumina": 2, "nextera": 1, "small_rna": 0}[name]))
    return identities, probe_counts, selected, content.hexdigest(), [digest.hexdigest() for digest in mate_content]


def command_output(command, label):
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, encoding="utf-8", errors="replace")
    require(result.returncode == 0,
            f"{label}: command failed with exit {result.returncode}: {shlex.join(command)}\n{result.stdout.strip()}")
    return result.stdout.strip()


def production_evidence(read1, read2):
    """Task logs are opportunistic evidence, not an extra public input contract."""
    resolved = [Path(path).resolve() for path in (read1, read2)]
    parents = list(dict.fromkeys(path.parent for path in resolved))
    files = []
    observed = []
    script_lines = []
    for parent in parents:
        for filename in (".command.log", ".command.out", ".command.err", ".command.sh", ".command.run"):
            path = parent / filename
            if not path.is_file():
                continue
            text = path.read_text(errors="replace")
            relevant = [line for line in text.splitlines()
                        if re.search(r"trim_galore|Trim Galore|trimgalore|adapter|selected|version|sha256:", line, re.I)]
            files.append({"path": str(path), "sha256": checksum(path), "relevant_lines": relevant})
            if filename in (".command.log", ".command.out", ".command.err"):
                observed.extend(relevant)
            else:
                script_lines.extend(relevant)
    return {
        "component": "trim_galore",
        "input_output_paths": {"rust_read1": str(resolved[0]), "rust_read2": str(resolved[1])},
        "task_files_available": bool(files),
        "task_files": files,
        "observed_version_lines": list(dict.fromkeys(line for line in observed if re.search(r"(?:version|Trim Galore).*2\.3\.0|trim_galore\s+2\.3\.0", line, re.I))),
        "observed_adapter_lines": list(dict.fromkeys(line for line in observed if re.search(r"adapter|selected", line, re.I))),
        "observed_invocation_lines": list(dict.fromkeys(line for line in observed if "--paired" in line and "trim_galore" in line)),
        "task_script_lines": list(dict.fromkeys(script_lines)),
        "availability_note": "Task evidence captured where mounted from the separate production-component parity rerun, not the original pipeline intermediates; integration must independently inspect rerun task logs and image digest when unavailable here.",
    }


def baseline(read1, read2, work_dir):
    executable = shutil.which("trim_galore_perl")
    require(executable is not None, "baseline runtime: trim_galore_perl is not installed")
    version = command_output([executable, "--version"], "Perl Trim Galore version")
    require(re.search(r"\b0\.6\.11\b", version), f"baseline runtime: Perl Trim Galore version is not 0.6.11: {version}")
    cutadapt = command_output(["cutadapt", "--version"], "Cutadapt version")
    require(cutadapt == "5.2", f"baseline runtime: Cutadapt version is not 5.2: {cutadapt}")
    output_dir = work_dir / "perl"
    output_dir.mkdir()
    command = [executable, "-q", "20", "--stringency", "2", "-e", "0.1", "--paired", "--gzip",
               "-o", str(output_dir), str(Path(read1).absolute()), str(Path(read2).absolute())]
    log_path = work_dir / "perl.log"
    with (
        log_path.open("w", buffering=1) as log,
        subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         text=True, bufsize=1, encoding="utf-8", errors="replace") as process,
    ):
        for line in process.stdout:
            print(line, end="", flush=True)
            log.write(line)
    require(process.returncode == 0, f"Perl baseline: exited {process.returncode}; inspect {log_path}")
    outputs = []
    for mate in (1, 2):
        matches = sorted(output_dir.glob(f"*_val_{mate}.fq.gz"))
        require(len(matches) == 1, f"Perl baseline: expected one validated R{mate} gzip output, got {len(matches)}")
        outputs.append(matches[0])
    log = log_path.read_text()
    selected = re.findall(r"^Adapter sequence: '([ACGT]+)'", log, re.M)
    require(selected and len(set(selected)) == 1, "Perl baseline: missing/ambiguous automatically selected adapter evidence")
    actual_adapters = re.findall(r"adapter sequence: '([ACGT]+)' from file .+? <<<", log)
    require(len(actual_adapters) == 2, "Perl baseline: missing per-mate Cutadapt invocation evidence")
    return outputs, {
        "tool": "official Perl Trim Galore",
        "version": "0.6.11",
        "version_output": version,
        "cutadapt_version": cutadapt,
        "perl_version_output": command_output(["perl", "-e", "printf \"%vd\\n\", $^V"], "Perl version"),
        "python_version": platform.python_version(),
        "source_url": PERL_SOURCE,
        "source_path": executable,
        "source_sha256": checksum(executable),
        "invocation": command,
        "shell_invocation": shlex.join(command),
        "adapter_selection": "independent baseline automatic selection; no preset supplied",
        "selected_adapter_r1": actual_adapters[0],
        "selected_adapter_r2": actual_adapters[1],
        "autodetected_adapter": selected[0],
        "adapter_log_lines": [line for line in log.splitlines() if re.search(r"adapter|Illumina\t|Nextera\t|smallRNA\t", line, re.I)],
        "log_path": str(log_path),
        "log_sha256": checksum(log_path),
    }


def compare_outputs(baseline_paths, rust_paths, identities):
    paired_digests = {name: hashlib.sha256() for name in ("perl", "rust")}
    mate_digests = {name: [hashlib.sha256(), hashlib.sha256()] for name in ("perl", "rust")}
    count = 0
    previous = {"perl": -1, "rust": -1}
    for position, (expected, actual) in enumerate(zip_longest(pairs(*baseline_paths, "Perl output"), pairs(*rust_paths, "Rust output")), 1):
        require(expected is not None and actual is not None, f"parity pair {position}: retained/missing/extra pair count mismatch")
        for tool, record in (("perl", expected), ("rust", actual)):
            require(record[0] in identities, f"{tool} output pair {position}: changed or extra identity {record[0]}")
            source_position = identities[record[0]]["position"]
            require(record[5:7] == identities[record[0]]["headers"],
                    f"{tool} output pair {position}: full header/comment changed for {record[0]}")
            require(source_position > previous[tool], f"{tool} output pair {position}: identity order differs from extracted inputs")
            previous[tool] = source_position
            record_digest(paired_digests[tool], record[:5])
            record_digest(mate_digests[tool][0], record[:3])
            record_digest(mate_digests[tool][1], (record[0], record[3], record[4]))
        require(expected[0] == actual[0], f"parity pair {position}: retained identity/order mismatch: {expected[0]} / {actual[0]}")
        for index, field in ((1, "R1 sequence"), (2, "R1 quality"), (3, "R2 sequence"), (4, "R2 quality")):
            require(expected[index] == actual[index], f"parity pair {position} ({expected[0]}): {field} mismatch")
        count += 1
    require(count > 0, "parity: no retained pairs")
    return count, {name: digest.hexdigest() for name, digest in paired_digests.items()}, {
        name: [digest.hexdigest() for digest in digests] for name, digests in mate_digests.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("id", "extracted_read1", "extracted_read2", "rust_read1", "rust_read2", "output", "work_dir"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    require(re.fullmatch(r"[^@\s]\S*", args.id), "sample ID: invalid whitespace or leading @")
    work_dir = Path(args.work_dir).absolute()
    work_dir.mkdir(parents=True, exist_ok=True)
    extracted = [args.extracted_read1, args.extracted_read2]
    rust = [args.rust_read1, args.rust_read2]
    identities, counts, selected, input_digest, input_mate_digests = inspect_inputs(*extracted)
    perl, perl_evidence = baseline(*extracted, work_dir)
    require(perl_evidence["selected_adapter_r1"] == PROBES[selected],
            "adapter selection: independent Perl baseline differs from the legacy three-probe decision")
    selected_r2 = "GATCGTCGGACT" if selected == "small_rna" else PROBES[selected]
    require(perl_evidence["selected_adapter_r2"] == selected_r2,
            "adapter selection: independent Perl baseline R2 adapter differs from the selected preset")
    retained, paired_digests, mate_digests = compare_outputs(perl, rust, identities)
    production = production_evidence(*rust)
    evidence = {
        "schema_version": 1,
        "sample_id": args.id,
        "passed": True,
        "complete_extracted_files": True,
        "comparison": "ordered (identity,R1 sequence,R1 quality,R2 sequence,R2 quality); identity=(first header token,terminal CB:Z: comment); full per-mate headers/comments must remain unchanged; gzip bytes and plus-line formatting ignored",
        "digest_encoding": "SHA-256 of newline-delimited compact ASCII JSON tuples",
        "counts": {"input_pairs": len(identities), "perl_retained_pairs": retained, "rust_retained_pairs": retained,
                   "discarded_pairs": len(identities) - retained},
        "mismatches": {"retained_identities": 0, "missing_identities": 0, "extra_identities": 0, "duplicates": 0,
                       "order": 0, "mate_names": 0, "changed_input_names_or_barcodes": 0, "full_headers_or_comments": 0,
                       "r1_sequence": 0, "r1_quality": 0, "r2_sequence": 0, "r2_quality": 0},
        "content_digests": {"input_pairs": input_digest, "input_mates": input_mate_digests,
                            "output_pairs": paired_digests, "output_mates": mate_digests},
        "files": {name: {"path": str(Path(path).absolute()), "sha256": checksum(path)}
                  for name, path in zip(("extracted_read1", "extracted_read2", "rust_read1", "rust_read2", "perl_read1", "perl_read2"), extracted + rust + perl)},
        "legacy_adapter_selection": {"contributing_records": min(len(identities), 999999), "probe_counts": counts,
                                     "selected_preset": "--" + selected, "adapter_r1": PROBES[selected], "adapter_r2": selected_r2,
                                     "basis": "recomputed from complete extracted R1 using the legacy contributing limit; not supplied to the independent baseline"},
        "baseline": perl_evidence,
        "production": production,
        "checker_invocation": [sys.executable] + sys.argv,
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
    print(f"PASS trim_galore parity ({args.id}): complete input={len(identities)}, retained={retained}, preset=--{selected}; zero ordered content/name mismatches")


if __name__ == "__main__":
    try:
        main()
    except (CheckError, OSError, ValueError, EOFError) as error:
        print(f"FAIL trim_galore parity: {error}", file=sys.stderr)
        sys.exit(1)
