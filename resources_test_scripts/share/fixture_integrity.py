#!/usr/bin/env python3
"""Validate portable fixture payloads and publish them without overwriting a set."""

import argparse
import contextlib
import csv
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import sys
import tempfile

CHROM = "chr14"
START = 23661200
END = 24661200
PREDICATE = 'rname=="chr14" && pos>=23661200 && pos<=24661200 && mpos>=23661200 && mpos<=24661200'
SAMPLE_IDS = ("H0O260120CPL0071F1S_001", "H0O260120CPL0071F1S_002")
REFERENCE_LENGTH = 107043718
REFERENCE_MD5 = "acbd9552c059d9b403e75ed26c1ce5bc"
REFERENCE_SOURCE = "/lustre/grp/xxllab/zhangzhen/FOODIE_Training_CP_20260305/resources/GRCh38/genome.fa"
CELL_LIST_SOURCES = {
    sample: Path("/lustre/grp/xxllab/fangj/03.project/02.snFOODIE/H0O260120CPL0071F1S/results/bam_cell_count/H0O260120CPL0071F1S") / sample / f"{sample}.bam.QC.cell.lst"
    for sample in SAMPLE_IDS
}
INTERVAL = {"chrom": CHROM, "start": START, "end": END,
            "coordinates": "1-based inclusive alignment and mate starts", "predicate": PREDICATE}
FIXTURE_SETS = ("common", "consensus_dedup", "raw_fastq_align")
REPO = Path(__file__).resolve().parents[2]


class FixtureError(ValueError):
    pass


def resources_root(repo, argument=None):
    return Path(os.path.abspath(Path(repo) / "resources_test" if argument is None else Path(argument).expanduser()))


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path, value):
    with open(path, "x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _probe_publication(stage):
    with tempfile.TemporaryDirectory(prefix=".publication-probe-", dir=stage) as temporary:
        source = Path(temporary) / "source"
        linked = Path(temporary) / "linked"
        with open(source, "xb") as handle:
            handle.write(b"fixture publication probe\n")
        try:
            os.link(source, linked)
            try:
                os.link(source, linked)
            except FileExistsError:
                pass
            else:
                raise FixtureError("filesystem did not refuse an existing payload link")
        except OSError as exc:
            raise FixtureError(
                f"{stage.parent}: fixture publication requires exclusive directory creation "
                f"and same-filesystem hard links: {exc}"
            ) from exc


def _publish_fixture(stage, destination, manifest):
    try:
        destination.mkdir()
    except FileExistsError as exc:
        raise FixtureError(f"refusing existing final destination: {destination}; explicitly remove it before rebuilding") from exc
    try:
        for name in sorted(manifest["files"]):
            target = destination / name
            target.parent.mkdir(parents=True, exist_ok=True)
            os.link(stage / name, target)
        # Readers must validate the set; only the marker makes it complete.
        os.link(stage / "completion.json", destination / "completion.json")
    except BaseException:
        print(f"Fixture publication failed or was interrupted at {destination}; validate before use. "
              "An incomplete final directory must be explicitly removed before rebuilding.", file=sys.stderr)
        raise


@contextlib.contextmanager
def staging_directory(destination, *, cleanup_on_error=True):
    destination = Path(os.path.abspath(destination))
    if os.path.lexists(destination):
        raise FixtureError(f"refusing existing final destination: {destination}; explicitly remove it before rebuilding")
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{destination.name}.incomplete-", dir=destination.parent))
    try:
        _probe_publication(stage)
        yield stage
        manifest = verify_fixture_set(destination.name, destination.parent, directory=stage)
        _publish_fixture(stage, destination, manifest)
        shutil.rmtree(stage)
    finally:
        if stage.exists():
            if cleanup_on_error:
                shutil.rmtree(stage)
            else:
                print(f"Retained private fixture stage: {stage}; confirm remote writers stopped before cleanup", file=sys.stderr)


def required_payloads(fixture_set):
    if fixture_set == "common":
        names = {"genome/chr14.fa"}
        for conversion in ("CT", "GA"):
            base = f"genome/Bisulfite_Genome/{conversion}_conversion"
            names.add(f"{base}/genome_mfa.{conversion}_conversion.fa")
            names.update(f"{base}/BS_{conversion}.{part}.bt2" for part in ("1", "2", "3", "4", "rev.1", "rev.2"))
        names.update(f"cell_lists/{sample}.cell_list.tsv" for sample in SAMPLE_IDS)
        return names
    if fixture_set == "consensus_dedup":
        return {"samples.csv", *(f"{sample}.bam" for sample in SAMPLE_IDS)}
    if fixture_set == "raw_fastq_align":
        return {"samples.csv", *(f"{sample}.{suffix}" for sample in SAMPLE_IDS
                                for suffix in ("R1.fastq.gz", "R2.fastq.gz", "identities.tsv"))}
    raise FixtureError(f"unsupported fixture set: {fixture_set!r}")


def _nonempty(path):
    if path.is_symlink() or not path.is_file() or path.stat().st_size == 0:
        raise FixtureError(f"missing, empty, or symlink payload: {path}")


def verify_reference(genome_dir):
    genome = Path(genome_dir)
    fasta = genome / "chr14.fa"
    _nonempty(fasta)
    digest = hashlib.md5()
    length = 0
    with open(fasta, "rb") as handle:
        header = handle.readline()
        if not header.startswith(b">chr14 ") and header.rstrip() != b">chr14":
            raise FixtureError(f"reference is not full chr14: {fasta}")
        for line in handle:
            sequence = line.strip()
            if sequence.startswith(b">"):
                raise FixtureError(f"unexpected additional reference sequence: {fasta}")
            digest.update(sequence.upper())
            length += len(sequence)
    if length != REFERENCE_LENGTH or digest.hexdigest() != REFERENCE_MD5:
        raise FixtureError(f"reference identity mismatch: {fasta}, length={length}, md5={digest.hexdigest()}")
    for name in required_payloads("common"):
        if name.startswith("genome/"):
            _nonempty(genome / name.removeprefix("genome/"))
    return {"name": CHROM, "length": length, "md5": digest.hexdigest(),
            "path": str(fasta.resolve()), "source": REFERENCE_SOURCE}


def _read_json(path):
    def unique_keys(pairs):
        value = {}
        for key, item in pairs:
            if key in value:
                raise FixtureError(f"duplicate JSON key {key!r}: {path}")
            value[key] = item
        return value
    try:
        with open(path, encoding="utf-8") as handle:
            value = json.load(handle, object_pairs_hook=unique_keys)
    except (OSError, ValueError) as exc:
        raise FixtureError(f"cannot read valid JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise FixtureError(f"expected a JSON object: {path}")
    return value


def _hash_map(value, expected, context):
    if not isinstance(value, dict) or set(value) != expected:
        raise FixtureError(f"{context}: required payload coverage mismatch")
    if any(not isinstance(digest, str) or re.fullmatch(r"[0-9a-f]{64}", digest) is None for digest in value.values()):
        raise FixtureError(f"{context}: invalid SHA-256 digest")


def _payloads(directory, expected):
    allowed_dirs = {str(parent) for name in expected for parent in Path(name).parents if str(parent) != "."}
    actual = set()
    for path in directory.rglob("*"):
        name = path.relative_to(directory).as_posix()
        if path.is_symlink():
            raise FixtureError(f"symlinks are not fixture payloads: {path}")
        if path.is_dir():
            if name not in allowed_dirs:
                raise FixtureError(f"unexpected owned directory: {path}")
        elif name != "completion.json":
            actual.add(name)
    if actual != expected:
        raise FixtureError(f"{directory}: owned payload mismatch; missing={sorted(expected - actual)}, extra={sorted(actual - expected)}")
    for name in expected:
        _nonempty(directory / name)


def _verify_csv(fixture_set, directory):
    fields = (["id", "input", "cell_list"] if fixture_set == "consensus_dedup"
              else ["id", "read1", "read2", "cell_list"])
    with open(directory / "samples.csv", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != fields:
            raise FixtureError(f"{directory}: samples.csv columns must be {fields}")
        rows = list(reader)
    if [row.get("id") for row in rows] != list(SAMPLE_IDS):
        raise FixtureError(f"{directory}: samples.csv must contain the two fixed IDs in order")
    for sample, row in zip(SAMPLE_IDS, rows):
        expected = {"id": sample, "cell_list": f"../common/cell_lists/{sample}.cell_list.tsv"}
        if fixture_set == "consensus_dedup":
            expected["input"] = f"{sample}.bam"
        else:
            expected.update(read1=f"{sample}.R1.fastq.gz", read2=f"{sample}.R2.fastq.gz")
        if row != expected:
            raise FixtureError(f"{directory}: invalid samples.csv paths or fields for {sample}")


def _missing_common(root):
    generate = shlex.join(["bash", str(REPO / "resources_test_scripts/cluster/common.sh"), "--resources-dir", str(root)])
    download = shlex.join(["bash", str(REPO / "resources_test_scripts/local/download.sh"), "common", "--resources-dir", str(root)])
    return f"common is missing under selected resources root {root}. On cluster login run: {generate}\nOr locally run: {download}"


def verify_fixture_set(fixture_set, resources_dir, *, directory=None):
    root = resources_root(REPO, resources_dir)
    selected = root / fixture_set if directory is None else Path(directory).absolute()
    expected = required_payloads(fixture_set)
    try:
        common = verify_fixture_set("common", root) if fixture_set != "common" else None
        if fixture_set == "common" and directory is None and not os.path.lexists(selected):
            raise FixtureError(_missing_common(root))
        if selected.is_symlink() or not selected.is_dir():
            raise FixtureError(f"fixture directory missing or symlink: {selected}")
        _nonempty(selected / "completion.json")
        manifest = _read_json(selected / "completion.json")
        fields = {"schema_version", "fixture_set", "completed", "created_at", "samples", "files", "dependencies"}
        if fixture_set == "consensus_dedup":
            fields.add("interval")
        if set(manifest) != fields:
            raise FixtureError(f"{selected}: invalid completion fields (legacy layouts must be regenerated)")
        if (type(manifest["schema_version"]) is not int or manifest["schema_version"] != 3
                or manifest["fixture_set"] != fixture_set or manifest["completed"] is not True
                or manifest["samples"] != list(SAMPLE_IDS)):
            raise FixtureError(f"{selected}: invalid completion schema/set/status/ordered samples")
        timestamp = manifest["created_at"]
        if not isinstance(timestamp, str):
            raise FixtureError(f"{selected}: created_at must be a UTC timestamp")
        try:
            parsed = datetime.datetime.fromisoformat(timestamp.replace("Z", "+00:00"))
            if parsed.utcoffset() != datetime.timedelta(0):
                raise ValueError("not UTC")
        except ValueError as exc:
            raise FixtureError(f"{selected}: invalid UTC created_at") from exc
        _hash_map(manifest["files"], expected, str(selected))
        _payloads(selected, expected)
        for name, digest in manifest["files"].items():
            if sha256_file(selected / name) != digest:
                raise FixtureError(f"fixture checksum mismatch: {selected / name}")
        if fixture_set == "common":
            if manifest["dependencies"] != {}:
                raise FixtureError(f"{selected}: common must have no dependencies")
            verify_reference(selected / "genome")
            for sample in SAMPLE_IDS:
                path = selected / f"cell_lists/{sample}.cell_list.tsv"
                with open(path, encoding="utf-8") as handle:
                    if not any(line.split("\t", 1)[0].strip() for line in handle):
                        raise FixtureError(f"empty listed-barcode content: {path}")
        else:
            dependencies = manifest["dependencies"]
            if not isinstance(dependencies, dict) or set(dependencies) != {"common"}:
                raise FixtureError(f"{selected}: expected exactly the common dependency")
            dependency = dependencies["common"]
            if not isinstance(dependency, dict) or set(dependency) != {"files"}:
                raise FixtureError(f"{selected}: invalid common dependency snapshot")
            _hash_map(dependency["files"], required_payloads("common"), f"{selected} common dependency")
            if dependency["files"] != common["files"]:
                raise FixtureError(f"{selected}: common payload changed; restore compatible common or explicitly remove/rebuild this dependent set")
            _verify_csv(fixture_set, selected)
            if fixture_set == "consensus_dedup":
                interval = manifest["interval"]
                if (not isinstance(interval, dict) or interval != INTERVAL
                        or type(interval.get("start")) is not int or type(interval.get("end")) is not int):
                    raise FixtureError(f"{selected}: invalid inclusive interval metadata")
        return manifest
    except (OSError, ValueError, KeyError, TypeError, EOFError) as exc:
        message = f"{selected}: {exc}"
        if fixture_set == "common" and os.path.lexists(root / "common"):
            message += f"\nCommon under {root} is invalid; restore compatible common, or explicitly remove/rebuild common and affected dependent sets. Do not rerun a generator over an existing destination."
        raise FixtureError(message) from exc


def common_cell_lists(resources_dir):
    root = resources_root(REPO, resources_dir)
    verify_fixture_set("common", root)
    return {sample: root / f"common/cell_lists/{sample}.cell_list.tsv" for sample in SAMPLE_IDS}


def write_completion(fixture_set, directory, *, common_manifest=None):
    directory = Path(directory)
    expected = required_payloads(fixture_set)
    _payloads(directory, expected)
    dependencies = {}
    if fixture_set != "common":
        if not isinstance(common_manifest, dict) or common_manifest.get("fixture_set") != "common":
            raise FixtureError("a validated common manifest is required")
        _hash_map(common_manifest.get("files"), required_payloads("common"), "common snapshot")
        dependencies = {"common": {"files": dict(common_manifest["files"])}}
    manifest = {"schema_version": 3, "fixture_set": fixture_set, "completed": True,
                "created_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                "samples": list(SAMPLE_IDS), "files": {name: sha256_file(directory / name) for name in sorted(expected)},
                "dependencies": dependencies}
    if fixture_set == "consensus_dedup":
        manifest["interval"] = dict(INTERVAL)
    write_json(directory / "completion.json", manifest)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    validate = sub.add_parser("validate")
    validate.add_argument("fixture_set", choices=FIXTURE_SETS)
    validate.add_argument("--resources-dir")
    args = parser.parse_args()
    try:
        root = resources_root(REPO, args.resources_dir)
        verify_fixture_set(args.fixture_set, root)
    except FixtureError as exc:
        parser.exit(1, f"error: {exc}\n")
    print(f"Validated {args.fixture_set}: {root / args.fixture_set}")


if __name__ == "__main__":
    main()
