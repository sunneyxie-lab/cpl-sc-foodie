# Test fixtures

## Ownership and roots

The three sets are independent transfer units:

| Set | Owned payloads | Prerequisite |
| --- | --- | --- |
| `common` | Full chr14 FASTA, both converted FASTAs, 12 Bowtie2 index parts, two unchanged original cell lists (17 files) | Cluster reference and lists |
| `consensus_dedup` | Two regional BAMs and `samples.csv` | Valid `common` |
| `raw_fastq_align` | Two paired gzip FASTQs, identity TSVs and `samples.csv` (7 files) | Valid `common` |

`cluster/common.sh` orchestrates shared preparation on login. `cluster/generate/`
contains compute workers and workflow-specific generators. `cluster/upload.sh`
uploads from login, `local/download.sh` downloads locally, and `share/` contains
validation, raw extraction and transfer helpers. No scientific workflow is defined
in this directory.

All shell entrypoints accept `--resources-dir <root>`. The default is this
checkout's absolute `resources_test/`, regardless of your current directory.
Relative explicit roots resolve against your initial directory. Quote paths with
spaces. Nextflow tests retain the separate `--resources_test <root>` spelling.

Generate or download `common` first. Nothing automatically generates or transfers
a prerequisite. Each workflow needs only `common` and its own fixture set. Both
CSVs reference `../common/cell_lists/<id>.cell_list.tsv`, relative to the CSV's
directory. Genomes come directly from `<root>/common/genome`, without a CSV column.
The sample order is `H0O260120CPL0071F1S_001`, then `H0O260120CPL0071F1S_002`.

Both regional selections use this exact samtools expression:

```text
rname=="chr14" && pos>=23661200 && pos<=24661200 && mpos>=23661200 && mpos<=24661200
```

These are 1-based inclusive alignment and mate starts. This isn't an overlap query
or an alignment-end constraint. The reference and cell lists aren't cropped.
Full chr14 has length `107043718` and normalized sequence MD5
`acbd9552c059d9b403e75ed26c1ce5bc`. Consensus generation reports distinct listed
barcode coverage and requires positive overlap, plus at least two BAM barcodes.
It no longer requires 90% input-list coverage. Scientific output checks remain
unchanged, including at least two output barcodes and the no-list cell-calling run.

## Prepare locally, then copy to the cluster

Use Viash 0.9.7, Nextflow 25.10.6 and Java 17 through 24. Build locally:

```bash
viash ns build --parallel --setup cachedbuild
```

The normal build includes the three checker modules used by local integration:
`target/_test/nextflow/test_workflows/{consensus_dedup_test,raw_fastq_align_test,trim_galore_parity_test}/main.nf`.
Their source scope stays `test`; don't push their images.

Copy the checkout's scripts, `src/` configs, `assets/`, the complete generated
`target/nextflow/` tree and any referenced local dependency resources to the cluster.
The preparation module needs its adjacent configuration and supporting files,
including `.config.vsh.yaml` and `nextflow_labels.config`. Don't copy only `main.nf`.
Don't install/build Viash on the cluster or download the full raw sources locally.
Registry image references don't require a Docker daemon on the cluster.

### Cluster access and SFTP

Open the bastion interactively:

```bash
ssh zhangzhen@172.18.18.11
```

Select `impetus集群登录节点 - gliese581` (`172.18.18.5`), then account
`zhangzhen`. Plain SSH remote-command execution through the bastion was rejected;
interactive menu selection worked.

Copy project files through SFTP:

```bash
sftp -P 9527 zhangzhen@172.18.18.4
```

An existing operator-authenticated SSH control connection can be reused by adding
`-o ControlPath=/path/to/control-socket`. The socket path is session-specific.

The project paths map as follows:

| Access | Project directory |
| --- | --- |
| SFTP | `/upload/pipelines/cpl-sc-foodie` |
| Cluster | `/lustre/sftp/xxllab/zhangzhen/upload/pipelines/cpl-sc-foodie` |

Transfer `src/`, `resources_test_scripts/`, `assets/`, the complete
`target/nextflow/` tree, referenced `target/dependencies/`, and root files
`_viash.yaml`, `nextflow.config`, `main.nf`, `CHANGELOG.md` and `.gitignore`.
Use `put -p` to preserve file modes and timestamps. Verify copied bytes by SFTP
readback and SHA-256 comparison, and check generator executable permissions.
Leave fixture data, credentials, Git metadata, runtime outputs and unrelated remote
files untouched. The fixture root must be shared between login and compute nodes
and pass the exclusive-directory/hard-link publication probe.

## Generate on the cluster

Choose an absolute shared root visible to compute nodes and login. Generators probe
exclusive directory creation and same-filesystem hard links before payload work.
The Lustre publication protocol doesn't require `RENAME_NOREPLACE`: after private
validation it exclusively creates the final directory, links each payload without
overwriting, then links `completion.json` last. The final directory is visible while
incomplete; consumers must validate it before use. There is no scratch-to-SFTP export.

Run common preparation **on login**, never by submitting the launcher with `sbatch`:

```bash
bash resources_test_scripts/cluster/common.sh --resources-dir /shared/fixtures
```

The login driver requires Python 3, Nextflow/Java, Singularity, `sbatch`, `scancel`,
`squeue`, readable original sources and the copied preparation module. It submits
one native shell extraction/list-copy job with `sbatch --parsable --wait` on
`cpu1,cpu2` (1 CPU, 4 GB, 8 hours). The logged job ID and successful wait status
precede preparation. The worker extracts the complete chr14 record and copies both
original lists byte-for-byte. It has no network, Nextflow or credential dependency.

Then one Nextflow controller stays on login and runs the existing Bismark preparation
module. `cluster/cluster.config` submits preparation and publish tasks to Slurm on
`cpu1,cpu2`, with an 8-hour task limit. Preparation labels request 10 CPUs/50 GB.
Nextflow pulls missing images on login with `ociAutoPull=false`; tasks use cached
native SIF images. The cluster's `singularity` command reports Apptainer 1.4.2.
Leave `ociMode` unset: explicit `false` makes Nextflow 25.10.6 emit `--no-oci`,
which this runtime rejects. The registry mirror applies to unqualified publish
images without rewriting the fully qualified GHCR tool image. No preload step is needed.

`NXF_SINGULARITY_CACHEDIR` defaults to `$HOME/.singularity/nxf_cache`. Its actual
location, the selected resource root and the run workspace must be shared with
compute nodes. Unset inherited `NXF_OFFLINE=true` if missing images need pulling.
Extraction and genome preparation run on compute; login only orchestrates, hashes
and publishes metadata. Hidden `.common.work-*` directories hold run work/logs,
outside the final payload.

After common succeeds, submit from the checkout (create `logs/` first):

```bash
mkdir -p logs
sbatch resources_test_scripts/cluster/generate/consensus_dedup.sh --resources-dir /shared/fixtures
sbatch resources_test_scripts/cluster/generate/raw_fastq_align.sh --resources-dir /shared/fixtures
```

Under `sbatch`, wrappers locate the checkout using `SLURM_SUBMIT_DIR`. Plain
allocation execution can use their installed paths. Consensus requests 16 CPUs,
32 GB and 8 hours, with two concurrent subsets and samtools >=1.12. Raw requests
16 CPUs, 32 GB and 24 hours, with Python 3, samtools >=1.12 and bqtools **0.6.1**.
Neither workflow generator needs Nextflow, Singularity, internet or cloud credentials.
Original BAM/FASTQ/list/reference paths remain the fixed cluster paths in the helpers.
`assets/barcodes.txt` remains repository-owned with its pinned SHA-256.

Raw extraction scans both original gzip streams through EOF, including after every
match. It preserves all four original lines, LF/CRLF endings and source order.
Real CBQ header grep supplies candidates; exact identities and original fields
are checked before materializing selected originals. Final gzip records are checked
against the original audit hashes before publication. No per-sample provenance is
distributed: later validation checks payload integrity, not the original-source
comparison. Private CBQ work, `counts.json` and `cbq.completed.json` remain under
`<root>/.raw_fastq_align_cbq`, or `RAW_FASTQ_ALIGN_CBQ_CACHE` when set, and are never
transferred.

## Validate and transfer one set

Python 3 is required on the workstation and Nextflow test controller. Validation
needs no raw-source tools or credentials:

```bash
python3 resources_test_scripts/share/fixture_integrity.py validate common --resources-dir /shared/fixtures
python3 resources_test_scripts/share/fixture_integrity.py validate consensus_dedup --resources-dir /shared/fixtures
python3 resources_test_scripts/share/fixture_integrity.py validate raw_fastq_align --resources-dir /shared/fixtures
```

Schema-3 completion certifies exact payload coverage, nonempty files and SHA-256
bytes. Workflow markers snapshot every common payload hash. Validation checks
common's actual files first, then compares the snapshot. Changing only common's
completion timestamp is compatible. Changing payload bytes invalidates existing
dependents. Restore compatible common or explicitly remove/rebuild affected sets;
transferring common doesn't validate or repair every dependent automatically.
Legacy markers/layouts and extra owned payloads are rejected.

### Convert existing raw fixtures

Stop fixture writers, transfers and readers. Before updating the checkout, validate
the old set with its old validator. Keep a backup of its `samples.csv`,
`completion.json` and both `.provenance.json` files outside the owned directory.
Then run from the updated checkout, using the same resources root:

```bash
python3 - /shared/fixtures <<'PY'
import csv
import io
import json
from pathlib import Path
import sys

sys.path.insert(0, "resources_test_scripts/share")
from fixture_integrity import SAMPLE_IDS, sha256_file, verify_fixture_set, write_completion

root = Path(sys.argv[1]).resolve()
directory = root / "raw_fastq_align"
common = verify_fixture_set("common", root)
marker = directory / "completion.json"
old = json.loads(marker.read_text())
if old["dependencies"]["common"]["files"] != common["files"]:
    raise SystemExit("common payload changed; restore the original dependency first")
preserved = {name: digest for name, digest in old["files"].items()
             if name.endswith((".fastq.gz", ".identities.tsv"))}
for name, digest in preserved.items():
    if sha256_file(directory / name) != digest:
        raise SystemExit(f"changed payload: {name}")
fields = ["id", "read1", "read2", "cell_list"]
with (directory / "samples.csv").open(newline="") as handle:
    rows = [{key: row[key] for key in fields} for row in csv.DictReader(handle)]
csv_text = io.StringIO()
writer = csv.DictWriter(csv_text, fieldnames=fields, lineterminator="\n")
writer.writeheader()
writer.writerows(rows)
marker.unlink()
(directory / "samples.csv").write_text(csv_text.getvalue())
for sample in SAMPLE_IDS:
    (directory / f"{sample}.provenance.json").unlink()
converted = write_completion("raw_fastq_align", directory, common_manifest=common)
verify_fixture_set("raw_fastq_align", root)
assert {name: converted["files"][name] for name in preserved} == preserved
print("Converted raw fixtures; FASTQ and identity TSV hashes unchanged.")
PY
```

On failure, restore the backed-up metadata with `completion.json` last, then
validate using the old checkout. Never regenerate a marker over changed payloads.
After success, remove the metadata backup. Convert the cluster copy and upload
`raw_fastq_align`; `rclone sync` removes the obsolete OSS provenance files. Download
that set to other workstations. Do not regenerate FASTQs or identity TSVs, change
`common`/`consensus_dedup`, or download an old-layout set with the updated validator.

Transfer commands require exactly one of `common`, `consensus_dedup` or
`raw_fastq_align`. `all`, multiple sets and obsolete options are rejected. The set
and `--resources-dir` may appear in either order. Install Python 3 and rclone.
Each invocation automatically sources its required readable credential file,
even when OSS variables already exist:

| Operation | Credential file |
| --- | --- |
| Upload on login | `~/.config/oss/xielab-readwrite.env` |
| Local download | `~/.config/oss/xielab-readonly.env` |

Each file must define `OSS_ENDPOINT`, `OSS_BUCKET`, `OSS_ACCESS_KEY_ID` and
`OSS_ACCESS_KEY_SECRET`. Optional `OSS_CA_CERT` must name a readable private CA.
Keep credentials outside the repository. TLS verification stays enabled, requests
use path-style S3 and bypass endpoint proxies, and rclone configuration stays in
the environment with bounded connection timeouts/retries.

Upload each named set separately from login, common first:

```bash
bash resources_test_scripts/cluster/upload.sh common --resources-dir /shared/fixtures
bash resources_test_scripts/cluster/upload.sh consensus_dedup --resources-dir /shared/fixtures
bash resources_test_scripts/cluster/upload.sh raw_fastq_align --resources-dir /shared/fixtures
```

Download locally, common first:

```bash
bash resources_test_scripts/local/download.sh common --resources-dir /local/fixtures
bash resources_test_scripts/local/download.sh consensus_dedup --resources-dir /local/fixtures
bash resources_test_scripts/local/download.sh raw_fastq_align --resources-dir /local/fixtures
```

Only `<root>/<set>/` and
`oss:$OSS_BUCKET/cpl-sc-foodie/resources_test/<set>/` are synchronized. Upload
validates before `rclone sync --progress` and replaces the selected remote prefix,
including deleting remote-only files. Download requires valid common for workflow
sets, syncs only the selected directory (deleting local-only files), then fully
validates it. Neither operation touches sibling sets or automatically transfers common.
Transfers update destinations in place, without transactional replacement or rollback.
A copied/old marker alone never proves successful transfer.

## Run scientific integration

```bash
src/workflows/consensus_dedup/integration_test.sh --resources_test /local/fixtures
src/workflows/raw_fastq_align/integration_test.sh --resources_test /local/fixtures
```

Both launchers retain Docker/no-publish execution and argument forwarding. Their
`test.nf` entrypoints invoke the shared Python validator synchronously before
constructing channels or submitting scientific tasks. Direct Nextflow invocation
cannot bypass this preflight. Validator helpers are declared as test resources at
`resources_test_scripts/share/` in both workflow configs.

Consensus expects three completed events: both supplied-list samples and the first
sample without a list. Raw expects both samples, complete Perl/Rust trimming parity
and unchanged first-pass BAM handoff. To prove fixture independence, run each with
its own private root containing only common and that workflow's set. Don't delete
another workflow's user data. Resume scientific Nextflow tasks with `-resume`.

## Failure and restart

- Every generator refuses an existing final path, including incomplete directories
  and dangling symlinks. Rebuilds require explicit operator removal. There's no reuse
  option or destructive reset. A restart always creates a fresh private sibling stage.
- Preflight/probe/extraction/preparation/validation failures never publish a final
  set. Common retains and reports its hidden stage/workspace even after cancellation.
  Its extraction job is cancelled by recorded ID, or by the unique run name scoped
  to the submitting user before ID capture. `squeue` confirms termination when
  available. Killing the local `sbatch --wait` client doesn't prove the remote job
  stopped. During preparation Nextflow cancels its own jobs. Confirm all writers
  stopped before explicitly cleaning retained paths. Scheduler/query failures don't
  authorize deletion.
- Raw/consensus remove only their own failed stage after local writers stop. They
  never remove existing fixtures, sources or another invocation's stage. Interrupted
  raw encoding can retain private CBQ work; only its encode-completion metadata
  certifies that archive, never the final regional set. Don't transfer partial gzip
  output or reuse it as a completed fixture.
- Timeout, SIGKILL or a lost launcher can leave hidden incomplete paths and surviving
  jobs. Inspect the logged extraction ID/name and Nextflow logs, then inspect/cancel
  surviving writers before cleanup. Errors remain in launcher/Slurm logs and exit
  status; there is no separate notification system.
- A validated private stage isn't a published fixture. Publication reserves the final
  directory exclusively and publishes its completion marker only after all payloads.
  A concurrent publisher loses without replacing the winner. Failure during linking
  leaves an incomplete final directory without a valid completion marker; after
  confirming writers stopped, explicitly remove that directory before rebuilding.
  The generator never removes it automatically. If interruption follows the last
  marker link, the final set may already be complete; validate it before deciding
  whether cleanup or rebuilding is needed. Never promote a failed stage by hand.
- Interrupted sync or failed post-download validation can leave the selected
  destination invalid. Don't use it. Rerun that named sync and validation. If common
  bytes changed, restore compatible common or rebuild dependents explicitly.

## Verification boundary

Local synthetic checks cover real samtools/bqtools extraction, full-reference
identity, schema/dependency corruption, publication races, worker isolation and
controller rejection. Scheduler/Nextflow stand-ins cover common transition rules;
a local rclone stand-in covers selected-set mirror/deletion semantics. Those harnesses
don't prove Slurm placement, login-side image pulling, actual cloud payload transfer
or scientific success on newly generated real fixtures.

Real acceptance still requires a supported shared cluster root, source access,
Slurm/Singularity on login, upload credentials, fresh generation, three selected
uploads/downloads, and both scientific integrations. Retain the fixed interval and
scientific assertions if smaller fixtures fail; report the observed failure rather
than weakening coverage. No speedup is claimed without measurement.
