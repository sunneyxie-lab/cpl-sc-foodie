#!/usr/bin/env python3
"""Generate and independently validate the cluster's paired consensus BAM subsets."""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import contextlib
import csv
import os
from pathlib import Path
import re
import shlex
import signal
import threading
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "share"))
from fixture_integrity import (
    CHROM,
    END,
    FixtureError,
    PREDICATE,
    SAMPLE_IDS,
    START,
    resources_root,
    staging_directory,
    verify_fixture_set,
    write_completion,
)

BAM_ROOT = Path("/lustre/grp/xxllab/fangj/03.project/02.snFOODIE/H0O260120CPL0071F1S/reRun/Align_bam")
BARCODE_SUFFIX = re.compile(r"_CB:Z:([^_]+)\Z")


def require_samtools():
    try:
        result = subprocess.run(
            ["samtools", "--version"], check=True, capture_output=True, text=True
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise FixtureError("stage=preflight: samtools is required on PATH") from exc
    match = re.match(r"samtools (\d+)\.(\d+)", result.stdout)
    if not match or tuple(map(int, match.groups())) < (1, 12):
        version = result.stdout.splitlines()[0] if result.stdout else "unknown version"
        raise FixtureError(f"stage=preflight: {version}; samtools >=1.12 is required for view -e")


@contextlib.contextmanager
def sam_records(bam):
    """Stream SAM records and reap the reader even when validation raises."""
    process = subprocess.Popen(
        ["samtools", "view", str(bam)], stdout=subprocess.PIPE, text=True
    )
    try:
        yield process.stdout
    except BaseException:
        process.kill()
        raise
    finally:
        process.stdout.close()
        status = process.wait()
    if status:
        raise FixtureError(f"samtools view failed while validating {bam} (exit {status})")


def generate_consensus_sample(source_bam, cell_list, destination, threads):
    """Subset explicit source paths, then validate actual output rather than the command."""
    source_bam, cell_list, destination = map(Path, (source_bam, cell_list, destination))
    sample_id = source_bam.stem
    if os.path.lexists(destination):
        raise FixtureError(f"sample={sample_id} stage=preflight: refusing existing destination: {destination}")
    if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
        raise FixtureError(f"sample={sample_id} stage=preflight: threads must be a positive integer")
    owned_output = False
    step = "preflight"
    try:
        for path in (source_bam, cell_list):
            if not path.is_file() or not os.access(path, os.R_OK):
                raise FixtureError(f"cannot read input {path}")
        with open(cell_list, encoding="utf-8") as handle:
            listed = {line.rstrip("\r\n").split("\t", 1)[0].strip() for line in handle}
        listed.discard("")
        if not listed:
            raise FixtureError(f"cell list has no listed barcodes: {cell_list}")

        step = "subset"
        print(f"{sample_id}: subsetting to {CHROM}:{START}-{END}", flush=True)
        # Exclusive creation also protects direct calls from racing over the same output.
        with open(destination, "xb") as output:
            owned_output = True
            subprocess.run(
                ["samtools", "view", "-h", "-b", "--no-PG", "-@", str(threads),
                 "-e", PREDICATE, str(source_bam)],
                stdout=output, check=True,
            )
        step = "quickcheck"
        subprocess.run(["samtools", "quickcheck", str(destination)], check=True)

        step = "validate-records"
        names = Counter()
        barcodes = set()
        record_count = 0
        with sam_records(destination) as records:
            for number, line in enumerate(records, 1):
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 11:
                    raise FixtureError(f"record {number}: expected at least 11 SAM fields")
                qname = fields[0]
                try:
                    position, mate_position = int(fields[3]), int(fields[7])
                except ValueError as exc:
                    raise FixtureError(f"record {number}: invalid alignment/mate start") from exc
                if not (fields[2] == CHROM and START <= position <= END
                        and START <= mate_position <= END):
                    raise FixtureError(f"record {number} ({qname}): violates inclusive predicate {PREDICATE}")
                suffix = BARCODE_SUFFIX.search(qname)
                if suffix is None:
                    raise FixtureError(f"record {number} ({qname}): no terminal _CB:Z:<barcode> read-name suffix")
                names[qname] += 1
                barcodes.add(suffix.group(1))
                record_count += 1
        if not record_count:
            raise FixtureError("subset is empty")
        invalid_pairs = []
        for qname, count in names.items():
            if count != 2:
                invalid_pairs.append(f"{qname} ({count})")
                if len(invalid_pairs) == 5:
                    break
        if invalid_pairs:
            raise FixtureError("read names not occurring exactly twice: " + ", ".join(invalid_pairs))
        if len(barcodes) < 2:
            raise FixtureError(f"{len(barcodes)} distinct BAM barcodes, expected at least 2")

        step = "cell-list-overlap"
        found = len(listed & barcodes)
        percentage = 100 * found / len(listed)
        print(f"{sample_id}: {record_count} records, {len(barcodes)} barcodes; "
              f"{found}/{len(listed)} cell list barcodes ({percentage:.2f}%) occur in the subset",
              flush=True)
        if found == 0:
            raise FixtureError("no listed barcode occurs in the subset")
        return {"records": record_count, "barcodes": len(barcodes),
                "listed": len(listed), "found": found}
    except BaseException as exc:
        # All local subprocesses have exited before removing this call's output.
        if owned_output:
            destination.unlink(missing_ok=True)
        if isinstance(exc, (OSError, UnicodeError, subprocess.SubprocessError, FixtureError)):
            raise FixtureError(f"sample={sample_id} stage={step}: {exc}") from exc
        raise


def samtools_thread_budget():
    allocated = os.environ.get("SLURM_CPUS_PER_TASK")
    if allocated is not None:
        try:
            cpus = int(allocated)
        except ValueError as exc:
            raise FixtureError("stage=preflight: SLURM_CPUS_PER_TASK must be a positive integer") from exc
        if cpus < 1:
            raise FixtureError("stage=preflight: SLURM_CPUS_PER_TASK must be a positive integer")
    else:
        cpus = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    # The two subsets run concurrently; -@ counts threads beyond samtools' main thread.
    return cpus // 2 - 1 if cpus > 3 else 1


@contextlib.contextmanager
def sample_workers():
    # Defer catchable interruption while writers run, including executor shutdown.
    # Raising during Thread.join could otherwise let staging cleanup race a writer.
    pending_signals = []
    previous_handlers = {}
    if threading.current_thread() is threading.main_thread():
        for signum in (signal.SIGINT, signal.SIGTERM):
            previous_handlers[signum] = signal.signal(
                signum, lambda signum, frame: pending_signals.append(signum)
            )
    try:
        with ThreadPoolExecutor(max_workers=2) as workers:
            yield workers
    finally:
        for signum, handler in previous_handlers.items():
            signal.signal(signum, handler)
    if pending_signals:
        raise KeyboardInterrupt(f"received signal {pending_signals[0]}")


def generate_fixture_set(repo, resources_dir):
    repo = Path(repo).resolve()
    resources_dir = resources_root(repo, resources_dir)
    destination = resources_dir / "consensus_dedup"
    if os.path.lexists(destination):
        raise FixtureError(f"refusing existing final destination: {destination}")
    common_manifest = verify_fixture_set("common", resources_dir)
    require_samtools()
    threads = samtools_thread_budget()
    sources = {sample_id: BAM_ROOT / f"{sample_id}.bam" for sample_id in SAMPLE_IDS}
    for sample_id, source in sources.items():
        if not source.is_file() or not os.access(source, os.R_OK):
            raise FixtureError(f"sample={sample_id} stage=preflight: cannot read source {source}; run on the lab cluster")
    with staging_directory(destination) as stage:
        # Executor shutdown waits for every worker, including after failure/interruption.
        # The staging context therefore never removes paths with active local writers.
        with sample_workers() as workers:
            futures = [workers.submit(
                generate_consensus_sample, sources[sample_id],
                resources_dir / "common" / "cell_lists" / f"{sample_id}.cell_list.tsv",
                stage / f"{sample_id}.bam", threads,
            ) for sample_id in SAMPLE_IDS]
            for future in futures:
                future.result()
        with open(stage / "samples.csv", "x", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(["id", "input", "cell_list"])
            for sample_id in SAMPLE_IDS:
                writer.writerow([sample_id, f"{sample_id}.bam",
                                 f"../common/cell_lists/{sample_id}.cell_list.tsv"])
        write_completion("consensus_dedup", stage, common_manifest=common_manifest)
        verify_fixture_set("consensus_dedup", resources_dir, directory=stage)
    return destination


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resources-dir", metavar="ROOT", help="fixture root (default: checkout/resources_test)")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[3]
    root = resources_root(repo, args.resources_dir)

    def interrupted(signum, frame):
        raise KeyboardInterrupt(f"received signal {signum}")

    signal.signal(signal.SIGTERM, interrupted)
    try:
        destination = generate_fixture_set(repo, root)
    except KeyboardInterrupt:
        print("error: consensus generation interrupted; local workers stopped before stage cleanup", file=sys.stderr)
        return 130
    except (FixtureError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    command = shlex.join(["bash", str(repo / "resources_test_scripts/cluster/upload.sh"),
                          "consensus_dedup", "--resources-dir", str(root)])
    print(f"Fixtures ready in {destination}; upload from the login node with:\n{command}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
