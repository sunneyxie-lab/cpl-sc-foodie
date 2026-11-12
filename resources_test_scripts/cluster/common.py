#!/usr/bin/env python3
"""Build shared fixtures with direct Slurm extraction and login-side Nextflow."""

import argparse
import contextlib
import getpass
import os
from pathlib import Path
import re
import selectors
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import uuid

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "resources_test_scripts/share"))
from fixture_integrity import (
    CELL_LIST_SOURCES,
    FixtureError,
    REFERENCE_SOURCE,
    SAMPLE_IDS,
    resources_root,
    sha256_file,
    staging_directory,
    verify_fixture_set,
    verify_reference,
    write_completion,
)

# Slurm limits execution to eight hours; also bound submission/queue waiting.
EXTRACTION_WAIT_SECONDS = 24 * 60 * 60
CANCELLATION_WAIT_SECONDS = 120
JOB_ID = re.compile(r"([0-9]+)(?:;([A-Za-z0-9_.-]+))?")


def log(message):
    print(message, file=sys.stderr, flush=True)


def readable_file(path):
    if not path.is_file() or not os.access(path, os.R_OK) or path.stat().st_size == 0:
        raise FixtureError(f"cannot read nonempty required file: {path}")


def preflight(repo, destination):
    if os.path.lexists(destination):
        raise FixtureError(f"refusing existing final destination: {destination}")
    if os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_JOBID"):
        raise FixtureError("common.sh must run on cluster login, not inside sbatch or a compute allocation; Nextflow pulls images on login")
    if sys.version_info.major != 3:
        raise FixtureError("Python 3 is required")
    diagnostics = {
        "nextflow": "nextflow not found on PATH. Install it with: curl -s https://get.nextflow.io | NXF_VER=25.10.6 bash (needs Java 17-24)",
        "java": "Java not found on PATH; Nextflow 25.10.6 needs Java 17-24; unpack a JDK 17 (e.g. Temurin) and set JAVA_HOME",
        "singularity": "singularity not found on PATH; ask the cluster admins",
    }
    for command in ("nextflow", "java", "singularity", "sbatch", "scancel", "squeue"):
        java_home = Path(os.environ["JAVA_HOME"]) / "bin/java" if os.environ.get("JAVA_HOME") else None
        available = os.access(java_home, os.X_OK) if command == "java" and java_home else shutil.which(command)
        if not available:
            raise FixtureError(diagnostics.get(command, f"{command} not found on PATH; run on cluster login with Slurm available"))
    try:
        result = subprocess.run(["nextflow", "-version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=60)
    except subprocess.TimeoutExpired as exc:
        raise FixtureError("nextflow -version timed out; check Nextflow 25.10.6 and Java 17-24") from exc
    if result.returncode:
        raise FixtureError("nextflow -version failed. Nextflow 25.10.6 needs Java 17-24; if the cluster's Java is older, unpack a JDK 17 (e.g. Temurin) and set JAVA_HOME\n" + result.stdout)

    module = repo / "target/nextflow/bismark/bismark_prepare"
    for name in ("main.nf", "nextflow.config", "nextflow_labels.config", ".config.vsh.yaml"):
        path = module / name
        try:
            readable_file(path)
        except FixtureError as exc:
            raise FixtureError(f"complete generated preparation directory required: {path}; build target/nextflow/bismark/bismark_prepare/ on the workstation and copy the entire directory to this checkout") from exc
    for path in module.rglob("*"):
        if path.is_file():
            readable_file(path)
        elif not path.is_dir():
            raise FixtureError(f"unreadable generated preparation support path: {path}; copy the entire workstation-built directory")
    for path in (repo / "resources_test_scripts/cluster/cluster.config", repo / "src/workflows/utils/labels.config", repo / "resources_test_scripts/cluster/generate/common_extract.sh"):
        readable_file(path)
    for path in (Path(REFERENCE_SOURCE), *(Path(CELL_LIST_SOURCES[sample]) for sample in SAMPLE_IDS)):
        readable_file(path)

    cache = Path(os.environ.get("NXF_SINGULARITY_CACHEDIR") or Path.home() / ".singularity/nxf_cache").expanduser().resolve()
    try:
        cache.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryFile(dir=cache) as probe:
            probe.write(b"cache-write-probe")
            probe.flush()
    except OSError as exc:
        raise FixtureError(f"Singularity image cache is not writable: {cache}; set NXF_SINGULARITY_CACHEDIR to a writable shared directory visible on compute") from exc
    os.environ["NXF_SINGULARITY_CACHEDIR"] = str(cache)
    log(f"Shared image cache: {cache} (must also be visible on compute nodes)")
    if os.environ.get("NXF_OFFLINE", "").lower() == "true":
        log("Warning: inherited NXF_OFFLINE=true disables required image pulls; unset NXF_OFFLINE on login if images are missing")
    return module


@contextlib.contextmanager
def ignore_interruptions():
    previous = {sig: signal.signal(sig, signal.SIG_IGN) for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)}
    try:
        yield
    finally:
        for sig, handler in previous.items():
            signal.signal(sig, handler)


def stop_local_process(process):
    if process.poll() is None:
        process.terminate()
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()


def cancel_extraction(job, name, user):
    # The random name is unique to this invocation, including before ID capture.
    cluster = ["--clusters", job[1]] if job and job[1] else []
    identity = [job[0]] if job else ["--user", user, "--name", name]
    query = ["--jobs", job[0]] if job else ["--name", name]
    log(f"Cancelling extraction {job[0] if job else name}")
    try:
        result = subprocess.run(["scancel", *cluster, *identity], capture_output=True, text=True, timeout=30)
        if result.returncode:
            log(f"Warning: scancel failed: {result.stderr.strip()}; checking squeue regardless")
        deadline = time.monotonic() + CANCELLATION_WAIT_SECONDS
        while True:
            result = subprocess.run(["squeue", *cluster, "--noheader", "--user", user, *query, "--format", "%i"], capture_output=True, text=True, timeout=30)
            if result.returncode:
                log(f"Warning: squeue failed; extraction termination is unconfirmed: {result.stderr.strip()}")
                return False
            if not result.stdout.strip():
                log("Extraction termination confirmed by squeue")
                return True
            if time.monotonic() >= deadline:
                log("Warning: extraction is still present in squeue; retain all shared paths")
                return False
            time.sleep(2)
    except (OSError, subprocess.TimeoutExpired) as exc:
        log(f"Warning: scheduler cancellation/query failed; extraction termination is unconfirmed: {exc}")
        return False


def run_extraction(repo, stage, workspace, *, timeout=EXTRACTION_WAIT_SECONDS):
    name = "fixture-common-" + uuid.uuid4().hex
    user = getpass.getuser()
    job = None
    process = None
    command = [
        "sbatch", "--parsable", "--wait", "--partition", "cpu1,cpu2",
        "--cpus-per-task", "1", "--mem", "4G", "--time", "8:00:00",
        "--job-name", name, "--output", str(workspace / "extract.%j.out"),
        "--error", str(workspace / "extract.%j.err"),
        str(repo / "resources_test_scripts/cluster/generate/common_extract.sh"),
        str(REFERENCE_SOURCE), str(workspace / "chr14.fa"),
    ]
    for sample in SAMPLE_IDS:
        command.extend([str(CELL_LIST_SOURCES[sample]), str(stage / "cell_lists" / f"{sample}.cell_list.tsv")])
    log(f"Submitting extraction {name}; logs: {workspace}")
    pending = b""

    def parse_line(line):
        nonlocal job
        text = line.decode(errors="replace").strip()
        match = JOB_ID.fullmatch(text)
        if match:
            parsed = match.groups()
            if job is not None and job != parsed:
                raise FixtureError("sbatch returned conflicting extraction job IDs")
            job = parsed
            log(f"Extraction Slurm job ID: {job[0]}" + (f";{job[1]}" if job[1] else ""))
            (workspace / "extraction-job.txt").write_text(text + "\n")
        elif text:
            log(text)

    try:
        # A separate session keeps terminal interrupts under this driver's control.
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, start_new_session=True)
        deadline = time.monotonic() + timeout
        with selectors.DefaultSelector() as selector, (workspace / "sbatch.stdout.log").open("wb") as stdout_log, (workspace / "sbatch.stderr.log").open("wb") as stderr_log:
            selector.register(process.stdout, selectors.EVENT_READ, stdout_log)
            selector.register(process.stderr, selectors.EVENT_READ, stderr_log)
            while selector.get_map() or process.poll() is None:
                if time.monotonic() >= deadline:
                    raise FixtureError("timed out waiting for extraction submission/queue/execution; preparation was not started")
                for key, _ in selector.select(timeout=0.5):
                    chunk = os.read(key.fileobj.fileno(), 65536)
                    if not chunk:
                        selector.unregister(key.fileobj)
                        continue
                    key.data.write(chunk)
                    key.data.flush()
                    if key.fileobj is process.stdout:
                        pending += chunk
                        while b"\n" in pending:
                            line, pending = pending.split(b"\n", 1)
                            parse_line(line)
                    else:
                        log(chunk.decode(errors="replace").rstrip())
            if pending:
                parse_line(pending)
        status = process.wait()
        if status:
            raise FixtureError(f"extraction sbatch --wait exited {status}; preparation was not started; inspect {workspace}")
        if job is None:
            raise FixtureError("sbatch --parsable --wait returned no valid job ID; preparation was not started")
        log(f"Extraction job {job[0]} completed successfully")
    except BaseException:
        with ignore_interruptions():
            # Stop the submit client first so a delayed submission cannot race a
            # successful empty name query. The remote job is cancelled separately.
            if process is not None:
                stop_local_process(process)
            cancel_extraction(job, name, user)
        raise
    finally:
        if process is not None:
            process.stdout.close()
            process.stderr.close()


def run_preparation(repo, module, stage, workspace):
    command = [
        "nextflow", "-log", str(workspace / "nextflow.log"), "run", str(module / "main.nf"),
        "-c", str(module / "nextflow.config"), "-profile", "singularity",
        "-c", str(repo / "resources_test_scripts/cluster/cluster.config"),
        "-work-dir", str(workspace / "work"), "--input", str(workspace / "chr14.fa"),
        "--output", "genome", "--publish_dir", str(stage / ".prepare"),
    ]
    with (workspace / "nextflow.console.log").open("wb") as console:
        process = subprocess.Popen(command, cwd=workspace, stdout=console, stderr=subprocess.STDOUT, start_new_session=True)
        log(f"Running preparation on login with Slurm tasks; log: {workspace / 'nextflow.console.log'}")
        try:
            status = process.wait()
        except BaseException:
            with ignore_interruptions():
                # Nextflow handles SIGINT by cancelling its own submitted jobs.
                if process.poll() is None:
                    os.killpg(process.pid, signal.SIGINT)
                    try:
                        process.wait(timeout=CANCELLATION_WAIT_SECONDS)
                    except subprocess.TimeoutExpired:
                        log("Warning: Nextflow cancellation is unconfirmed; retaining all shared paths and logs")
            raise
    if status:
        advice = " Unset NXF_OFFLINE on login so Nextflow can pull missing images." if os.environ.get("NXF_OFFLINE", "").lower() == "true" else ""
        raise FixtureError(f"Nextflow preparation failed with exit {status}; inspect {workspace / 'nextflow.console.log'}.{advice}")
    prepared = stage / ".prepare/genome"
    if not prepared.is_dir() or prepared.is_symlink():
        raise FixtureError(f"successful Nextflow run did not publish a genome directory: {prepared}")
    prepared.rename(stage / "genome")
    # Standalone generated wrappers also publish state metadata, not payload.
    shutil.rmtree(stage / ".prepare")


def generate_common(repo, root):
    destination = root / "common"
    module = preflight(repo, destination)
    stage = None
    workspace = None
    try:
        with staging_directory(destination, cleanup_on_error=False) as stage:
            workspace = Path(tempfile.mkdtemp(prefix=".common.work-", dir=root))
            log(f"Private stage: {stage}; private workspace: {workspace}")
            (stage / "cell_lists").mkdir()
            run_extraction(repo, stage, workspace)
            for sample in SAMPLE_IDS:
                copied = stage / "cell_lists" / f"{sample}.cell_list.tsv"
                readable_file(copied)
                if sha256_file(copied) != sha256_file(CELL_LIST_SOURCES[sample]):
                    raise FixtureError(f"copied cell list differs from original: {sample}")
            readable_file(workspace / "chr14.fa")
            run_preparation(repo, module, stage, workspace)
            verify_reference(stage / "genome")
            write_completion("common", stage)
            verify_fixture_set("common", root, directory=stage)
    except BaseException:
        if stage is not None:
            log(f"Failed common stage retained for inspection and explicit cleanup: {stage}")
        if workspace is not None:
            log(f"Failed common workspace retained (remote writers may still be active): {workspace}")
        raise
    # Publication happens on context exit; all successful writers have finished.
    shutil.rmtree(workspace)
    return destination


def main():
    parser = argparse.ArgumentParser(description="Build common fixtures on cluster login; Slurm runs all sequence extraction/preparation")
    parser.add_argument("--resources-dir", metavar="ROOT")
    args = parser.parse_args()

    def interrupted(signum, frame):
        raise InterruptedError(f"interrupted by {signal.Signals(signum).name}")

    for sig in (signal.SIGTERM, signal.SIGHUP):
        signal.signal(sig, interrupted)
    try:
        destination = generate_common(REPO, resources_root(REPO, args.resources_dir))
    except (FixtureError, OSError, KeyboardInterrupt) as exc:
        log(f"error: {exc or 'interrupted'}")
        return 1
    log(f"Completed common fixtures: {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
